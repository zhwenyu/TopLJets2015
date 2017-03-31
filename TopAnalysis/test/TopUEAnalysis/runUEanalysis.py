#!/usr/bin/env/python

import glob
import sys
import os
import optparse
import ROOT
import numpy as np
import array as array
import pickle

from UEEventCounter import *
from UEAnalysisHandler import *
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist

#GLOBAL VARIABLES TO DEFINE THE ANALYSIS

BASEQUANTILES=[ 100*x/50 for x in xrange(0,51) ]

"""
determines the resolutions needed to define the migration matrices
"""
def determineSliceResolutions(opt):

    #build the chain
    t=ROOT.TChain('tue')
    for f in opt.input.split(','): 
        if len(f):
            t.AddFile(f)

    #loop the available events and fill resolution arrays for events passing the selection cuts
    varVals={}
    for var in VARS.keys() : varVals[var]=[[],[]]
    ue=UEEventCounter(ptthreshold=[float(x) for x in opt.ptThr.split(',')])
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%1000==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()

        #require a pure event selected at reco and gen levels
        passSel=(t.passSel&0x1)
        gen_passSel=t.gen_passSel
        if not passSel or not gen_passSel: continue

        #count particles in the event
        ue.count(t,debug=False)
        #ue.show()
        #raw_input()

        #fill resolution arrays
        for var in VARS.keys():
            isAngle=VARS[var][4]
            try:
                val=getattr(t,'gen_%s'%var)
                deltaVal=getattr(t,var)[0]-val
            except:
                val=getattr(ue,'gen_%s'%var) 
                deltaVal=getattr(ue,'rec_%s'%var) - val

            if isAngle : 
                val      = val*180./ROOT.TMath.Pi()
                deltaVal = ROOT.TVector2.Phi_mpi_pi(deltaVal)*180./ROOT.TMath.Pi()
                val      = ROOT.TMath.Abs(val)
            varVals[var][0].append( val )
            varVals[var][1].append( deltaVal )

    #prepare canvas for the resolution plots
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
    c=ROOT.TCanvas('c','c',550,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.12)
    c.SetBottomMargin(0.1)

    #determine resolution for 10 quantiles
    varResolutions={}
    for var in varVals:

        #use quantiles and extremes to define the 2D resolution histogram
        genvarQ=[]
        if var=='nj': 
            genvarQ = [0,1,2,5]
        else : 
            genvarQ = np.percentile( np.array(varVals[var][0]), BASEQUANTILES )
        dvarQ = np.percentile( np.array(varVals[var][1]), [2.5,97.5])
        h2d=ROOT.TH2F(var,
                      ';%s (gen. level);#Delta(reco-gen);'%VARS[var][0],
                      len(genvarQ)-1,
                      array.array('d',genvarQ),50,dvarQ[0],dvarQ[1])
        h2d.GetZaxis().SetNdivisions(5)
        for i in xrange(0,len(varVals[var][0])):
            h2d.Fill(varVals[var][0][i],varVals[var][1][i])
        
        #fit a gausian to determine resolution per quantile
        resGr=ROOT.TGraphErrors()
        resGr.SetName(var+'_resol')
        resGr.SetMarkerStyle(20)
        for xbin in xrange(1,h2d.GetNbinsX()+1):
            
            tmp=h2d.ProjectionY('tmp',xbin,xbin)

            total=tmp.Integral(0,tmp.GetNbinsX()+1)
            if total==0 : continue

            #require minimum stats for resolution fit
            if total>10 :
                
                mean,sigma=tmp.GetMean(),tmp.GetRMS()
                if not var in ['C','D','aplanarity','sphericity']:
                    tmp.Fit('gaus','MQ+')
                    gaus=tmp.GetFunction('gaus')
                    mean=gaus.GetParameter(1)
                    sigma=gaus.GetParameter(2)
                npts=resGr.GetN()
                resGr.SetPoint(npts,h2d.GetXaxis().GetBinCenter(xbin),mean)
                resGr.SetPointError(npts,0,sigma)

            #normalize entries in quantile
            for ybin in xrange(1,h2d.GetNbinsY()+1):
                val=h2d.GetBinContent(xbin,ybin)
                err=h2d.GetBinError(xbin,ybin)
                h2d.SetBinContent(xbin,ybin,100.*val/total)
                h2d.SetBinError(xbin,ybin,100.*err/total)
            tmp.Delete()

        #show plot
        c.Clear()
        h2d.Draw('colz')
        h2d.GetYaxis().SetTitleOffset(1.1)
        resGr.Draw('p')
        tex=ROOT.TLatex()        
        tex.SetTextFont(42)
        tex.SetTextSize(0.035)
        tex.SetNDC()
        tex.DrawLatex(0.1,0.96,'#bf{CMS} #it{simulation preliminary} %s'%VARS[var][0])
        tex.DrawLatex(0.8,0.96,'#sqrt{s}=13 TeV')
        logx=True if ('chavg' in var or 'chflux' in var or 'ptttbar' in var or 'sumpt' in var or 'ptpos' in var or 'ptll' in var) else False
        c.SetLogx(logx)
        c.Modified()
        c.Update()
        for ext in ['png','pdf']: c.SaveAs('%s/%s_resol.%s'%(opt.out,var,ext))

        varResolutions[var]=(h2d,resGr)


    #compute correlations at generator level
    ROOT.gStyle.SetPaintTextFormat("4.0f")
    hcorr=ROOT.TH2F('slicecorr',';Variable; Variable; Correlation (%)',len(varVals),0,len(varVals),len(varVals),0,len(varVals))
    hcorr.GetZaxis().SetNdivisions(5)
    xbin=0
    for var1 in varVals:
        xbin+=1
        hcorr.GetXaxis().SetBinLabel(xbin,VARS[var1][0])
        ybin=0
        for var2 in varVals:
            rho=np.corrcoef(varVals[var1][0],varVals[var2][0])[0][1]
            ybin+=1
            hcorr.GetYaxis().SetBinLabel(ybin,VARS[var2][0])
            hcorr.SetBinContent(xbin,ybin,rho*100)
    c.Clear()
    hcorr.Draw('colztext')
    hcorr.GetYaxis().SetTitleOffset(1.25)
    tex=ROOT.TLatex()        
    tex.SetTextFont(42)
    tex.SetTextSize(0.035)
    tex.SetNDC()
    tex.DrawLatex(0.1,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.DrawLatex(0.8,0.96,'#sqrt{s}=13 TeV')
    c.Modified()
    c.Update()
    for ext in ['png','pdf']: c.SaveAs('%s/slicecorr.%s'%(opt.out,ext))

    #all done, save to pickle file
    with open(os.path.join(opt.out,'sliceresolutions.pck'), 'w') as cachefile:
        pickle.dump(varResolutions, cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(hcorr,          cachefile, pickle.HIGHEST_PROTOCOL)


"""
defines the variables to slice the UE measurement and saves their reco quantiles
"""
def defineAnalysisBinning(opt):

    #read the resolutions
    varResolutions=None
    with open(os.path.join(opt.out,'sliceresolutions.pck'),'r') as cachefile:
        varResolutions = pickle.load(cachefile)

    #
    # VARIABLES
    # readout resolution curves and determine the bins for the observables
    #
    print 'Defining axes for the slice variables'
    varAxes={}
    for var in varResolutions:

        #special case for jet multiplicity
        if var=='nj':
            genBin=[0,1,2,3]
            recBin=[0,1,2,3]
        else:
             
            #get resolution map and quantiles       
            genBin=[0]
            resolGr=varResolutions[var][1]
            nSigmaForBins=2.0
            lastAcceptResol=resolGr.GetErrorY(0)
            for n in xrange(1,resolGr.GetN()-1):

                #center value, bias and resolution (don't correct bias)
                xgen_i,delta_i=ROOT.Double(0),ROOT.Double(0)        
                resolGr.GetPoint(n,xgen_i,delta_i)
                dx=xgen_i-genBin[-1]
                if dx<0 : continue
                if dx<nSigmaForBins*lastAcceptResol : continue
                genBin.append( xgen_i )
                lastAcceptResol=resolGr.GetErrorY(n-1)
            genBin.append( genBin[-1]+lastAcceptResol*2 )

            
            recBin=[0]
            for i in xrange(1,len(genBin)):
                dBin=genBin[i]-genBin[i-1]
                if i>1 : recBin.append( genBin[i-1]+0.25*dBin )
                recBin.append( genBin[i-1]+0.75*dBin )
            recBin.append( genBin[-1]+0.25*(genBin[-1]+genBin[-2]) )


            #special case for angular variables: override previous definition
            if VARS[var][4]:
                nbins=10
                if var=='phittbar': nbins=4
                delta=180./float(nbins)
                genBin=[i*delta     for i in xrange(0,nbins+1)]
                recBin=[i*delta*0.5 for i in xrange(0,2*nbins+1)]

        #save binning in histos
        varAxes[(var,False)] = ROOT.TAxis(len(genBin)-1,array.array('d',genBin))
        varAxes[(var,False)].SetName('%s_genSlices'%var)
        varAxes[(var,True)]  = ROOT.TAxis(len(recBin)-1,array.array('d',recBin))
        varAxes[(var,True)].SetName('%s_recSlices'%var)
        print '%30s'%('Bin definition for '+var)
        for level in [False,True]:
            nbins=varAxes[(var,level)].GetNbins()
            print '%30s'%'','[%3.1f,%3.1f] / %d @ %s'%(varAxes[(var,level)].GetBinLowEdge(1),
                                                       varAxes[(var,level)].GetBinUpEdge(nbins),
                                                       nbins,
                                                       'rec' if level else 'gen')
            print '%30s'%'',(recBin if level else genBin)
        print ''

    #
    # DEFINE MIGRATION MATRICES, GEN/REC LEVEL HISTOS
    #
    histos={}    
    print 'Saving gen/rec level histos and migration matrix templates'
    print 'Observables:',OBSVARS
    print 'Event axes:',EVAXES
    print 'Slice vars:',SLICEVARS

    #inclusive
    for var in OBSVARS:
        for a in EVAXES+['inc']:
            nbins={}
            for level in [False,True]:

                name='%s_%s%s'%('rec' if level else 'gen',a,var)            
                nbins[level]=varAxes[(var,level)].GetNbins()
                if a != 'inc' : nbins[level]*=3
                histos[ (var,a,level) ] = ROOT.TH1F(name,name,nbins[level],0,nbins[level])
                histos[ (var,a,level) ].SetDirectory(0)
            name='m_%s%s'%(a,var)
            histos[ (var,a) ] = ROOT.TH2F(name,name,nbins[False],0,nbins[False],nbins[True],0,nbins[True])
            histos[ (var,a) ].SetDirectory(0)

    #sliced
    for var in SLICEVARS:
        
        nslicebins={
            False:varAxes[(var,False)].GetNbins(),
            True:varAxes[(var,True)].GetNbins()
            }

        for obs in OBSVARS:

            nbins={}
            for level in [False,True]:
                nbinsObs=varAxes[(obs,level)].GetNbins()
                nbins[level]=nslicebins[level]*nbinsObs

                name='%s_%s_%s'%('rec' if level else 'gen',obs,var)
                histos[ (obs,level,var) ] = ROOT.TH1F(name,name,nbins[level],0,nbins[level])
                histos[ (obs,level,var) ].SetDirectory(0)

            name='m_%s_%s'%(obs,var)
            histos[ (obs,var) ] = ROOT.TH2F(name,name,nbins[False],0,nbins[False],nbins[True],0,nbins[True])
            histos[ (obs,var) ].SetDirectory(0)

    #all done, save to pickle file
    with open(os.path.join(opt.out,'analysiscfg.pck'), 'w') as cachefile:
        pickle.dump(varAxes,     cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(histos,      cachefile, pickle.HIGHEST_PROTOCOL)


"""
loops over a set of files with common name to fill the migration matrices
"""        
def runUEAnalysis(inF,outF,wgtIdx,varIdx,cfgDir):

    print '[runAnalysis] %s -> %s wgtIdx=%d,varIdx=%d'%(inF,outF,wgtIdx,varIdx)
    
    #configure from pickle file
    ueHandler=UEAnalysisHandler(os.path.join(cfgDir,'analysiscfg.pck'))

    #loop over the tree to fill histos
    ue=UEEventCounter(EVAXES)
    t=ROOT.TChain('tue')
    t.AddFile(inF)
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%100==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()

        #count particles
        ue.count(t,varIdx=varIdx)

        #selection flags
        gen_passSel        = t.gen_passSel
        passSel=((t.passSel>>varIdx) & 0x1)


        #fill histos (if sliceVar is None, non-sliced histos are filled)
        for sliceVar in [None]+SLICEVARS:
            
            #check if sliceVar is non null and configure values to use
            sliceVarVals=None
            try:
                if sliceVar!='chmult':
                    sliceVarVals=(sliceVar, getattr(t,'gen_%s'%sliceVar), getattr(t,sliceVar)[varIdx] )
                else:
                    sliceVarVals=('chmult', ue.gen_chmult, ue.rec_chmult)
            except:
                pass

            #loop over UE observables
            for obs in OBSVARS:
                ueHandler.fillInclusive(sliceVarVals=sliceVarVals,
                                        obs=obs,
                                        ue=ue,
                                        weight=t.weight[wgtIdx] if passSel else 0.,
                                        gen_passSel=gen_passSel,
                                        passSel=passSel)

                #loop over axes defining away/towards/transverse regions
                #for a in EVAXES:
                #        ueHandler.fillDifferential(sliceVarVals=sliceVarVals,
                #                                   obs=obs,
                #                                   a=a,
                #                                   ue=ue,
                #                                   weight=t.weight[wgtIdx],
                #                                   gen_passSel=gen_passSel,
                #                                   passSel=passSel)


    #save histos to ROOT file
    fOut=ROOT.TFile.Open(outF,'RECREATE')
    for h in ueHandler.histos:
        ueHandler.histos[h].SetDirectory(fOut)
        ueHandler.histos[h].Write()
    fOut.Close()

"""
Wrapper for when the analysis is run in parallel
"""
def runUEAnalysisPacked(args):
    try:
        fileName,outDir,varIdx,wgtIdx,cfgDir=args
        runUEAnalysis(fileName,outDir,varIdx,wgtIdx,cfgDir)
    except : 
        print 50*'<'
        print "  Problem with", name, "continuing without"
        print 50*'<'
        return False


def main():

    ROOT.gROOT.SetBatch(True)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',    dest='input',  help='input',                       default='MC13TeV_TTJets_dilpowheg_0.root',   type='string')
    parser.add_option(      '--ptThr', dest='ptThr',  help='ptThreshold gen,reco',        default='1.0,0.9',   type='string')
    parser.add_option('-s', '--step',  dest='step',   help='step',                        default=1,   type=int)
    parser.add_option('-w', '--wgt',   dest='wgtIdx', help='weight index to use',         default=0,   type=int)
    parser.add_option('-v', '--var',   dest='varIdx', help='calib index to use',          default=0,   type=int)
    parser.add_option('-j', '--jobs',  dest='jobs',   help='jobs to run in parallel',     default=1,   type=int)
    parser.add_option('-o', '--out',   dest='out',    help='output',                      default='./UEanalysis',   type='string')
    parser.add_option(      '--only',  dest='only',   help='csv list of tags to process', default='',  type='string')
    parser.add_option('-q', '--queue', dest='queue',  help='Batch queue to use [default: %default]', default='local')
    (opt, args) = parser.parse_args()

    onlyList=opt.only.split('v')
    os.system('mkdir -p %s'%opt.out)

    if opt.step==0:
        determineSliceResolutions(opt)

    if opt.step==1:
        defineAnalysisBinning(opt)

    if opt.step==2:

        #prepare output
        outDir=opt.out+'/analysis_%d_%d/Chunks'%(opt.wgtIdx,opt.varIdx)
        os.system('mkdir -p %s'%outDir)

        #create the tasklist
        file_list=[]
        if os.path.isdir(opt.input):
            for file_path in os.listdir(opt.input):
                if file_path.endswith('.root'):
                    file_list.append(os.path.join(opt.input,file_path))
        elif opt.input.startswith('/store/'):
            file_list = getEOSlslist(opt.input)
        elif '.root' in opt.input:
            file_list.append(opt.input)


        tasklist=[]
        for filename in file_list:
            baseFileName=os.path.basename(filename)
            tag,ext=os.path.splitext(baseFileName)
            if len(onlyList)>0:
                processThis=False
                for filtTag in onlyList:
                    if filtTag in tag:
                        processThis=True
                if not processThis : continue
            tasklist.append((filename,'%s/%s'%(outDir,baseFileName),opt.wgtIdx,opt.varIdx,opt.out))

        #run jobs locally
        if opt.queue=='local':
            if opt.jobs>1:
                print ' Running %d jobs in %d threads' % (len(tasklist), opt.jobs)
                import multiprocessing as MP
                pool = MP.Pool(opt.jobs)
                pool.map(runUEAnalysisPacked,tasklist)
            else:
                for fileName,outfile,wgtIdx,varIdx,cfgDir in tasklist:
                    runUEAnalysis(fileName,outfile,wgtIdx,varIdx,cfgDir)
        #submit jobs
        else:
            print ' Running %d jobs to %s'%(len(tasklist),opt.queue)
            cmsswBase=os.environ['CMSSW_BASE']
            for fileName,_,wgtIdx,varIdx,cfgDir in tasklist:
                localRun='python %s/src/TopLJets2015/TopAnalysis/test/TopUEAnalysis/runUEanalysis.py -i %s -o %s -q local -s 2 -w %d -v %d'%(cmsswBase,fileName,cfgDir,wgtIdx,varIdx)
                cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localRun)
                os.system(cmd)


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
