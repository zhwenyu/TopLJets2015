import ROOT
import os
import sys
import argparse
import pickle
from TopLJets2015.TopAnalysis.roofitTools import showFitResult,shushRooFit

#sigma=1pb distributed accross crossing angles 
#NB this does not sum to 1 as we don't use all crossing angles in the analysis
SIGNALXSECS={120:0.269,130:0.273,140:0.143,150:0.293}
VALIDLHCXANGLES=[120,130,140,150]

def getCsiAcceptanceCuts(url,mass):

    """reads the pickle file with the linear parameterization of the csi acceptances from simulation
    computes the acceptance at a given mass for each crossing angle"""

    if not url : return {}
    with open(url,'r') as f:
        csiaccParam=pickle.load(f)

    csiacc={}
    for xangle in VALIDLHCXANGLES:
        csiacc[xangle]=[]
        for rp in [23,123]:
            a,b=csiaccParam[(xangle,rp)]
            csiacc[xangle].append( (a+b*mass/1000.,0.18) )
    return csiacc

        
def addPDFwithParameterUncertainties(w,pdfName,addParamTag,model_params,setParamConst=False,nuisMax=3):

    """ adds a parametrization to affect parameters by their uncertainties and creates a new PDF using those """

    cmd='EDIT::{0}_{1}({0},'.format(pdfName,addParamTag)

    #change from fixed parameter to one affected by a nuisance to parameterize intrinsic uncertainty
    #x -> x + sigma*nuis
    iter = model_params.createIterator()
    var = iter.Next()
    while var :
        name=var.GetName()
        val=var.getVal()
        unc=var.getError()        
        paramName='{0}_{1}'.format(name,addParamTag)
        if setParamConst:
            w.factory('{0}[{1}]'.format(paramName,val))
            cmd+='{0}={1},'.format(name,paramName)
        else:
            w.factory("expr::{0}_param('max({1}+{2}*@0,1e-6)',{{{0}_nuis[0,-{3},{3}]}})".format(paramName,val,unc,nuisMax) )            
            cmd+='{0}={1}_param,'.format(name,paramName)
        var = iter.Next()

    cmd=cmd[0:-1]+')'
    w.factory(cmd)

    return pdfName+'_'+addParamTag

def parametrizeBackground(w,args,opt):

    """paraterizes the signal and adds to the workspace"""    

    totalBkg={}

    #base background PDF to be specified by category
    w.factory('RooCBShape::bkg_decay(mmiss,bkg_mean_decay[1000,2000],bkg_sigma_decay[100,2000],bkg_alpha_decay[3,5],bkg_n_decay[4,6])')
    w.factory('RooCBShape::bkg_turnon(mmiss,bkg_mean_turnon[600,1000],bkg_sigma_turnon[100,300],bkg_alpha_turnon[3,5],bkg_n_turnon[4,6])')
    w.factory("SUM::bkg(bkg_decay_frac[0.2,0.6]*bkg_decay,bkg_turnon)")

    #import signal events
    data=ROOT.TChain('data')
    for f in [os.path.join(opt.input,x) for x in os.listdir(opt.input) if 'Data13TeV' in x]:
        if 'MuonEG' in f : continue
        data.AddFile(f)

    #define final preselection cuts
    cuts=opt.presel
    if opt.csiacc:
        csiCuts ='csi1>%f && csi1<%f && '%opt.csiacc[opt.xangle][0]
        csiCuts+='csi2>%f && csi2<%f'%opt.csiacc[opt.xangle][1]
        cuts=csiCuts if len(cuts)==0 else '{0} && {1}'.format(cuts,csiCuts)

    #replace the RP information with the one from the opt.mix
    print '\t replacing reconstructed values from the RPs with',opt.mix,'to parametrize background'
    orig_ds=ROOT.RooDataSet("bkg_data_orig","bkg_data_orig",data,args,opt.presel)
    ds=ROOT.RooDataSet("bkg_data","bkg_data",args)
    for i in xrange(orig_ds.numEntries()):
        iargs=orig_ds.get(i)
        for x in ['csi1','csi2','mpp','mmiss']:
            iargs.setRealValue(x, iargs.getRealValue(opt.mix+x))
        ds.add(iargs)
    ds=ds.reduce(cuts)

    #loop over categories and fit background shape
    for icat in range(len(opt.categs)):
        
        #apply category cuts
        categCut=opt.categs[icat]
        if len(categCut)==0: 
            orig_rds=orig_ds
            rds=ds            
        else:
            orig_rds=orig_ds.reduce(ROOT.RooFit.Cut(categCut))
            rds=ds.reduce(ROOT.RooFit.Cut(categCut))

        #import dataset
        getattr(w,'import')(orig_rds if opt.unblind else rds, 
                            ROOT.RooFit.Rename('data_obs_%s_a%d_cat%d'%(opt.chTag,opt.xangle,icat)) )

        if opt.unblind:
            print '*'*50
            print 'Observed data has been imported *unblinded* for cat=',icat
            print '*'*50


        totalBkg[icat]=rds.sumEntries()
        
        w.pdf('bkg').fitTo(rds)
        model_params=w.pdf('bkg').getParameters(rds)
        pdfName=addPDFwithParameterUncertainties(w,'bkg','%s_a%d_cat%d'%(opt.chTag,opt.xangle,icat),model_params)
        
        #show fit result
        showFitResult(w=w,
                      varDesc=('mmiss','Missing mass [GeV]'),
                      xran=(0,2500),
                      data=rds,
                      pdfCompList=[(pdfName, '*turnon*', ROOT.kGray, 2),
                                   (pdfName, None,     ROOT.kBlue, 1)],
                      extraText='%s, #alpha=%d, %s'%(opt.chTitle,opt.xangle,categCut),
                      outName='%s/bkgpdf_%s_a%d_cat%d'%(opt.output,opt.chTag,opt.xangle,icat))

    print '\t total background:',totalBkg
    return totalBkg


def parametrizeSignal(w,args,opt):

    """paraterizes the signal and adds to the workspace"""    

    #number of signal events expected
    sigExp={}

    #base PDF to be specified by category
    
    #simple Crystal-Ball
    #w.factory("RooCBShape::sig_core(mmiss, sig_mean_core[%3.0f,%3.0f],sig_sigma_core[20,200], sig_alpha_core[1,5], sig_n_core[0,10])"%(opt.mass-100,opt.mass+100))

    #resonant core as Crystal-Ball(resonant) * Gaussian(resol)
    #w.factory("RooBreitWigner::sig_res_core(mmiss, sig_mean_core[%3.0f,%3.0f],sig_width_core[0,200])"%(opt.mass-100,opt.mass+100))        
    w.factory("RooCBShape::sig_res_core(mmiss, sig_res_mean_core[%3.0f,%3.0f],sig_res_sigma_core[0,200], sig_res_alpha_core[0.1,10], sig_res_n_core[0,10])"%(opt.mass-100,opt.mass+100))
    w.factory("RooGaussian::sig_resol_core(mmiss, sig_resol_mean_core[0],sig_resol_sigma_core[0,100])")
    w.factory("FCONV::sig_core(mmiss,sig_res_core,sig_resol_core)") ;

    #simple combinatorial shape
    w.factory("RooCBShape:sig_comb(mmiss,sig_mean_comb[400,1500],sig_sigma_comb[100,400], sig_alpha_comb[0,10], sig_n_comb[0,10])")

    #combinatorial based on the background shape
    #w.factory('RooCBShape::sig_decay(mmiss,sig_mean_decay[1000,2000],sig_sigma_decay[100,2000],sig_alpha_decay[1,5],sig_n_decay[4,6])')
    #w.factory('RooCBShape::sig_turnon(mmiss,sig_mean_turnon[600,1000],sig_sigma_turnon[100,300],sig_alpha_turnon[1,5],sig_n_turnon[4,6])')
    #w.factory("SUM::sig_comb(sig_decay_frac[0.0,1.0]*sig_decay,sig_turnon)")
        
    w.factory("SUM::sig(sig_core_frac[0.6,1.0]*sig_core,sig_comb)")

    #define final preselection cuts
    cuts=opt.presel
    if opt.csiacc:
        csiCuts ='csi1>%f && csi1<%f && '%opt.csiacc[opt.xangle][0]
        csiCuts+='csi2>%f && csi2<%f'%opt.csiacc[opt.xangle][1]
        cuts=csiCuts if len(cuts)==0 else '{0} && {1}'.format(cuts,csiCuts)

    #import signal events
    fIn=ROOT.TFile.Open(os.path.join(opt.input,opt.sig))
    data=fIn.Get('data')
    wgtds=ROOT.RooDataSet("wgtsig_data","wgtsig_data",data,args,cuts,'wgt')        
    ds=ROOT.RooDataSet("sig_data","sig_data",args,ROOT.RooFit.Import(data),ROOT.RooFit.Cut(cuts))
    fIn.Close()

    #loop over categories 
    for icat in range(len(opt.categs)):
            
        nWgt=wgtds.sumEntries()
        categCut=opt.categs[icat]
        if len(categCut)==0:
            nWgt=wgtds.sumEntries()
            rds=ds
        else:
            wgtrds=wgtds.reduce(ROOT.RooFit.Cut(categCut))
            nWgt=wgtrds.sumEntries()
            rds=ds.reduce(ROOT.RooFit.Cut(categCut))

        #set initial values for the combinatorial component 
        #from the ones fit to the background
        #let only average and width be adjusted for the combinatorial component
        #i.e. use exponential tails and from background fits
        #comb_params=w.pdf('sig_comb').getParameters(rds)
        #iter=comb_params.createIterator()
        #param = iter.Next()
        #while param :
        #    name=param.GetName()
        #    bkgName=name.replace('sig_','bkg_')+'_%s_a%d_cat%d_param'%(opt.chTag,opt.xangle,icat)
        #    w.var(name).setVal(w.function(bkgName).getVal())
        #    param = iter.Next()

        #run the fit
        w.pdf('sig').fitTo(rds)
        sigExp[icat]=nWgt*SIGNALXSECS[opt.xangle]*opt.lumi*opt.mu
        model_params=w.pdf('sig').getParameters(rds)

        pdfName=addPDFwithParameterUncertainties(w,'sig','%s_a%d_cat%d'%(opt.chTag,opt.xangle,icat),model_params)

        #show fit result
        showFitResult(w=w,
                      varDesc=('mmiss','Missing mass [GeV]'),
                      xran=(0,2500),
                      data=rds,
                      pdfCompList=[(pdfName, '*comb*', ROOT.kGray, 2),
                                   (pdfName, None,     ROOT.kBlue, 1)],
                      extraText='%s, #alpha=%d, %s'%(opt.chTitle,opt.xangle,categCut),
                      outName='%s/sigpdf_%s_a%d_cat%d'%(opt.output,opt.chTag,opt.xangle,icat))
        
    print '\t signal expectations:',sigExp
    return sigExp


def writeDataCard(w,opt,sigExp,bkgExp):

    """writes the datacard and the workspace"""

    #write workspace to ROOT file
    wsROOT='workspace_%s_a%d.root'%(opt.chTag,opt.xangle)
    w.writeToFile(os.path.join(opt.output,wsROOT))

    #create a card per category
    dcList=[]
    for icat in range(len(opt.categs)):
        cat='%s_a%d_cat%d'%(opt.chTag,opt.xangle,icat)
        dcTxt='%s/shapes-parametric.datacard_%s.dat'%(opt.output,cat)
        dcList.append(dcTxt)
        with open(dcTxt,'w') as dc:
            dc.write('#\n')
            dc.write('# datacard was automatically generated with generateWorkspace.py\n')
            dc.write('# the options passed are printed below\n')
            dc.write('# %s\n'%opt)
            dc.write('#\n')
            dc.write('imax *\n')
            dc.write('jmax *\n')
            dc.write('kmax *\n')
            dc.write('-'*50+'\n')
            dc.write('shapes sig      {0} {1} w:$PROCESS_{0}\n'.format(cat,wsROOT))
            dc.write('shapes bkg      {0} {1} w:$PROCESS_{0}\n'.format(cat,wsROOT))
            dc.write('shapes data_obs {0} {1} w:data_obs_{0}\n'.format(cat,wsROOT))
            dc.write('-'*50+'\n')
            dc.write('bin %s\n'%cat)
            dc.write('observation -1\n')
            dc.write('-'*50+'\n')
            dc.write('%15s %15s %15s\n'%('bin',cat,cat))
            dc.write('%15s %15s %15s\n'%('process','sig','bkg'))
            dc.write('%15s %15s %15s\n'%('process','0', '1'))
            dc.write('%15s %15s %15s\n'%('rate','%3.2f'%sigExp[icat], '%3.2f'%bkgExp[icat]))
            dc.write('-'*50+'\n')
            
            #float the background normalization as well as the signal
            dc.write('mu_bkg{0} rateParam {0} bkg 1\n'.format(cat))

            args=w.allVars()
            it = args.createIterator()
            while it:
                try:                    
                    obj=it.Next()                    
                    vname = obj.GetName()
                    if 'nuis' in vname:
                        addToDataCard=True
                        if ('bkg' in vname or 'sig' in vname) and not str(cat) in vname:
                            addToDataCard=False            
                        if addToDataCard:
                            dc.write('%40s param 0 1\n'%vname)                    
                except:
                    break
            del it

    #write a local wrapper to run combine
    with open(os.path.join(opt.output,'runScan.sh'),'w') as script:
        script.write('#!/bin/bash\n')
        allDC=' '.join(['cat%d=%s'%(i,os.path.basename(dcList[i])) for i in range(len(dcList))])
        script.write('cd {0}\n'.format(opt.output))
        script.write('combineCards.py {0} > combined_card.dat\n'.format(allDC))
        script.write('text2workspace.py combined_card.dat\n')
        script.write('combine combined_card.dat.root -M AsymptoticLimits -t -1\n')
        script.write('combine combined_card.dat.root -M Significance -t -1 --expectSignal=1\n')
        script.write('cd -\n'.format(opt.output))
        
    print '\t workspace is available @',wsROOT
    print '\t generated the following datacards',dcList
    print '\t runScan.sh can be found in the same directory'


def main(args):

    parser = argparse.ArgumentParser(description='usage: %prog [options]')
    parser.add_argument('-i', '--input',
                        dest='input',   
                        default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/ab05162/analysis/',
                        help='input directory with the files [default: %default]')
    parser.add_argument('--xangle',
                        dest='xangle',
                        default=150,
                        type=int,
                        help='crossing angle [%default]')
    parser.add_argument('--sig',
                        dest='sig',
                        default='MC13TeV_ppxz_m800_x150.root',
                        help='signal point [%default]')
    parser.add_argument('--presel',
                        dest='presel', 
                        default='xangle==150 && cat==169 && l1pt>30 && l2pt>20 && bosonpt>50',
                        help='preselection [default: %default]')
    parser.add_argument('--csiacc',
                        dest='csiacc', 
                        default='test/analysis/pps/signal_resol_acc.pck',
                        help='parametrization of the csi acceptance cuts and mass resolution [default: %default]')
    parser.add_argument('--categs',
                        dest='categs',
                        default='nvtx<20,nvtx>=20',
                        help='Sub-categories [default: %default]')
    parser.add_argument('--mix',
                        dest='mix',
                        default='mix',
                        help='Mixing values to use [default: %default]')
    parser.add_argument('--lumi',
                        dest='lumi',
                        default=37500.,
                        type=float,
                        help='integrated luminosity [default: %default]')
    parser.add_argument('--mu',
                        dest='mu',
                        default=1.0,
                        type=float,
                        help='signal strength [default: %default]')
    parser.add_argument('-o', '--output',
                        dest='output', 
                        default='analysis/stat',
                        help='Output directory [default: %default]')
    parser.add_argument('--unblind',
                        dest='unblind', 
                        default=False,
                        action='store_true',
                        help='Use non-mixed data in the final fit [default: %default]')
    opt=parser.parse_args(args)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    shushRooFit()

    #parse channels from pre-selection string
    import re
    ch_dict={'169':'zmm','121':'zee','22':'g'}
    chtit_dict={'169':'Z#rightarrow#mu#mu','121':'Z#rightarrowee','22':'#gamma'}
    regex=re.compile('cat==(\d+)') 
    allCh=regex.findall(opt.presel)
    setattr(opt,'chTag','_'.join([ch_dict[x] for x in allCh]))
    setattr(opt,'chTitle',','.join([chtit_dict[x] for x in allCh]))

    #configuration
    if not opt.xangle in VALIDLHCXANGLES : 
        print 'Crossing angle',opt.xangle,'is not valid a valid one'
        return
    opt.categs=opt.categs.split(',')
    if len(opt.categs)==0 : opt.categs=[]
    opt.mass=int(opt.sig.split('_')[2].replace('m','')) #hardcoded...
    if opt.csiacc: opt.csiacc=getCsiAcceptanceCuts(opt.csiacc,opt.mass)

    #start the output
    os.system('mkdir -p %s'%opt.output)

    print '[generateWorkspace]'
    print '\t signal from',opt.sig,'mass=',opt.mass
    print '\t will apply the following preselection:',opt.presel
    print '\t channel tag=',opt.chTag
    if len(opt.categs) : 
        print '\t will categorize in:',opt.categs
    if opt.csiacc:
        print '\t will apply csi acceptance cuts',opt.csiacc

    #start the workspace
    w=ROOT.RooWorkspace('w')
    for pfix in ['','mix','mixem']:
        w.factory(pfix+'csi1[0,1]')
        w.factory(pfix+'csi2[0,1]')
        w.factory(pfix+'mpp[0,5000]')
        w.factory(pfix+'mmiss[0,5000]')
    w.factory('l1pt[20,6500]')
    w.factory('l2pt[20,6500]')
    w.factory('bosonpt[20,6500]')
    w.factory('xangle[100,200]')
    w.factory('nvtx[0,100]')
    w.factory('rho[0,1000]')
    w.factory('PFMultSumHF[0,100000]')
    w.factory('PFHtSumHF[0,100000]')
    w.factory('PFPzSumHF[0,100000]')
    w.factory('cat[0,200]')
    w.factory('wgt[-1.0e9,1.0e9]')
    varSet=w.allVars()

    #parametrize the background
    print '\t starting background parametrization'
    bkgExp=parametrizeBackground(w, varSet , opt)

    #parametrize the signal
    print '\t starting signal parametrization'
    sigExp=parametrizeSignal(w, varSet, opt)

    #write summary in datacards
    print '\t writing datacard'
    writeDataCard(w,opt,sigExp,bkgExp)
                          
    #move everything to output
    print '\t all done, output can be found in',opt.output

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))



