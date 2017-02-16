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

#var name, var title, use to slice phase space, use as event axis, is angle
VARS={
    'ptttbar'  : ('p_{T}(t#bar{t})',  True,  False, False),
    'phittbar' : ('#phi(t#bar{t})',   True,  True,  True),
    'ptpos'    : ('p_{T}(l^{+})',     True,  False, False),
    'phipos'   : ('#phi(l^{+})',      True,  True,  True),
    'ptll'     : ('p_{T}(l,l)',       True,  False, False),
    'phill'    : ('#phi(ll)',         True,  True,  True),
    'sumpt'    : ('#Sigma p_{T}(l)',  True,  False, False),
    'mll'      : ('M(l,l)',           True,  False, False),
    'dphill'   : ('#Delta#phi(l,l)',  True,  False, True),
    'nj'       : ('N(jets)',          True,  False, False),
    'chmult'   : ('N(ch)',            True,  False, False),
    'chflux'   : ('#Sigma p_{T}(ch)', False, False, False),
    'chavgpt'  : ('#bar{p}_{T}(ch)',  False, False, False),
    }


"""
determines the resolutions needed to define the migration matrices
"""
def determineSliceResolutions(opt):

    #build the chain
    t=ROOT.TChain('tue')
    for f in opt.input.split(','): t.AddFile(f)

    #loop the available events and fill resolution arrays for events passing the selection cuts
    varVals={}
    for var in VARS.keys() : varVals[var]=[[],[]]
    ue=UEEventCounter(ptthreshold=opt.ptThr)
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
        ue.count(t) #,debug=True)

        #fill resolution arrays
        for var in VARS.keys():
            isAngle=VARS[var][3]
            try:
                val=getattr(t,'gen_%s'%var)
                deltaVal=getattr(t,var)[0]-val
            except:
                val=getattr(ue,'gen_%s'%var) 
                deltaVal=getattr(ue,'rec_%s'%var) - val

            if isAngle : 
                val  = val*180./ROOT.TMath.Pi()
                deltaVal = ROOT.TVector2.Phi_mpi_pi(deltaVal)*180./ROOT.TMath.Pi()
                val  = ROOT.TMath.Abs(val)
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
                      ';%s (gen. level);Resolution;'%VARS[var][0],
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
                tmp.Fit('gaus','MQ+')
                gaus=tmp.GetFunction('gaus')
                npts=resGr.GetN()
                resGr.SetPoint(npts,h2d.GetXaxis().GetBinCenter(xbin),gaus.GetParameter(1))
                resGr.SetPointError(npts,0,gaus.GetParameter(2))

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
        hcorr.GetXaxis().SetBinLabel(xbin,VARS[var][0])
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
    # SLICE VARIABLES
    # readout resolution curves and determine the bins for the observables
    #
    print 'Defining axes for the slice variables'
    slicingAxes={}
    for var in varResolutions:

        #special case for jet multiplicity
        if var=='nj':
            genBin=[0,1,2,3]
            recBin=[0,1,2,3]
        else:
             
            #get resolution map and quantiles       
            inigenBin=[]           
            genResol=[]
            resolGr=varResolutions[var][1]
            x,y=ROOT.Double(0),ROOT.Double(0)        
            for n in xrange(0,resolGr.GetN()):
                resolGr.GetPoint(n,x,y)
                ey=resolGr.GetErrorY(n)
                inigenBin.append(float(x))
                genResol.append(float(ey))

            #use max. resolution as limiting factor
            maxRes=np.amax(genResol)
            #special case for angular variables
            if var=='phill' or var=='dphill' or var=='phittbar' or var=='phipos':             
                inigenBin=[i*18. for i in xrange(0,11)]
                maxRes=18.

            #define gen bin merging quantiles if needed according to the resolution
            i=0
            genBin=[ ROOT.TMath.Min(0,inigenBin[0]) ]
            if var=='ptttbar': genBin.append( inigenBin[0] )
            while True:

                if i>len(inigenBin)-2 : break
                
                #next bin defined is saturated by the resolution
                xi=inigenBin[i]
                xj=genBin[-1]
                xii=inigenBin[i+1]                
                dx=xii-xi
                if xii<=xj :
                    i+=1 
                    continue

                if dx>maxRes:
                    genBin.append(xii)
                else:
                    genBin.append(ROOT.TMath.Min(xj+maxRes,inigenBin[-1]))

            #define rec bins from genBins/2
            recBin=[]
            for i in xrange(1,len(genBin)):
                dx=(genBin[i]-genBin[i-1])*0.5
                recBin.append( genBin[i-1] )
                recBin.append( genBin[i-1]+dx)
            recBin.append( genBin[-1] )
        
        #save binning in histos
        slicingAxes[(var,False)] = ROOT.TAxis(len(genBin)-1,array.array('d',genBin))
        slicingAxes[(var,False)].SetName('%s_genSlices'%var)
        slicingAxes[(var,True)]  = ROOT.TAxis(len(recBin)-1,array.array('d',recBin))
        slicingAxes[(var,True)].SetName('%s_recSlices'%var)


    #
    # OBSERVABLES
    # use quantiles to determine the binning for the observables : nch, sumpt and avgpt
    #
    print 'Defining axes for the observables'
    obsVals={}
    for obs in OBSQUANTILES: 
        for level in [False,True]:
            obsVals[(obs,level)]={'inc':[[]]}
            for a in EVENTAXES: 
                obsVals[(obs,level)][a]=[[],[],[]]
    ue=UEEventCounter(EVENTAXES)
    t=ROOT.TChain('tue')
    for f in opt.input.split(','): t.AddFile(f)
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%100==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()

        #require a pure event selected at reco and gen levels
        passSel=(t.passSel&0x1)
        gen_passSel=t.gen_passSel
        if not passSel or not gen_passSel: continue
        ue.count(t)

        obsVals[('chmult',False)]['inc'][0].append(ue.gen_chmult)
        obsVals[('chflux',False)]['inc'][0].append(ue.gen_chflux)
        obsVals[('chavgpt',False)]['inc'][0].append(ue.gen_chavgpt) 
        obsVals[('chmult',True)]['inc'][0].append(ue.rec_chmult)
        obsVals[('chflux',True)]['inc'][0].append(ue.rec_chflux)
        obsVals[('chavgpt',True)]['inc'][0].append(ue.rec_chavgpt)
        for a in EVENTAXES:
            for k in xrange(0,3):
                obsVals[('chmult',False)][a][k].append( ue.gen_chmult_wrtTo[a][k] )
                obsVals[('chflux',False)][a][k].append( ue.gen_chflux_wrtTo[a][k] )
                obsVals[('chavgpt',False)][a][k].append( ue.gen_chavgpt_wrtTo[a][k] )
                obsVals[('chmult',True)][a][k].append( ue.rec_chmult_incWrtTo[a][k] )
                obsVals[('chflux',True)][a][k].append( ue.rec_chflux_incWrtTo[a][k] )
                obsVals[('chavgpt',True)][a][k].append( ue.rec_chavgpt_incWrtTo[a][k] )

    #determine quantiles and save as binnings
    obsAxes={}
    for key in obsVals:
        obs,level=key
        for a in obsVals[key]:
            for i in xrange(0,len(obsVals[key][a])):
                levelStr='rec' if level else 'gen'
                varQ = np.percentile( np.array( obsVals[key][a][i] ), OBSQUANTILES[obs])
                varQ[0]=0.                
                if level:
                    recBin=[]
                    for k in xrange(1,len(varQ)):
                        dx=(varQ[k]-varQ[k-1])*0.5
                        recBin.append( varQ[k-1] )
                        recBin.append( varQ[k-1]+dx)
                    recBin.append( varQ[-1] )
                    obsAxes[ (obs,a,i,level) ] = ROOT.TAxis(len(recBin)-1,array.array('d',recBin))
                else:
                    obsAxes[ (obs,a,i,level) ] = ROOT.TAxis(len(varQ)-1,array.array('d',varQ))
                obsAxes[ (obs,a,i,level) ].SetName('%s_%s%s_%d'%(levelStr,obs,a,i))


    #
    # DEFINE MIGRATION MATRICES, GEN/REC LEVEL HISTOS
    #
    print 'Saving gen/rec level histos and migration matrix templates'
    histos={}    
    for key in obsVals:
        if key[1] : continue
        obs=key[0]
        for a in obsVals[key]:

            nbins={}
            for level in [False,True]:
                levelStr='rec' if level else 'gen' 
                name='%s_%s%s'%(levelStr,a,obs)            
                nbins[level]=0
                newKey=(obs,level)
                for i in xrange(0,len(obsVals[newKey][a])): nbins[level] += obsAxes[ (obs,a,i,level) ].GetNbins()
                histos[ (obs,a,level) ] = ROOT.TH1F(name,name,nbins[level],0,nbins[level])
                histos[ (obs,a,level) ].SetDirectory(0)

            name='m_%s%s'%(a,obs)
            histos[ (obs,a) ] = ROOT.TH2F(name,name,nbins[False],0,nbins[False],nbins[True],0,nbins[True])
            histos[ (obs,a) ].SetDirectory(0)

    #sliced
    for var,isRec in slicingAxes:
        if isRec: continue

        nslicebins={
            False:slicingAxes[(var,False)].GetNbins(),
            True:slicingAxes[(var,True)].GetNbins()
            }

        for key in obsVals:
            if key[1]: continue
            obs=key[0]
            for a in obsVals[key]:

                nbins={}
                for level in [False,True]:

                    newKey=(obs,level)

                    nbinsObs=0
                    for i in xrange(0,len(obsVals[newKey][a])): nbinsObs += obsAxes[ (obs,a,i,level) ].GetNbins()
                    nbins[level]=nslicebins[level]*nbinsObs

                    levelStr='rec' if level else 'gen'
                    name='%s_%s%s_%s'%(levelStr,a,obs,var)
                    histos[ (obs,a,level,var) ] = ROOT.TH1F(name,name,nbins[level],0,nbins[level])
                    histos[ (obs,a,level,var) ].SetDirectory(0)

                name='m_%s%s_%s'%(a,obs,var)
                histos[ (obs,a,var) ] = ROOT.TH2F(name,name,nbins[False],0,nbins[False],nbins[True],0,nbins[True])
                histos[ (obs,a,var) ].SetDirectory(0)

    #all done, save to pickle file
    with open(os.path.join(opt.out,'analysiscfg.pck'), 'w') as cachefile:
        pickle.dump(obsAxes,     cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(slicingAxes, cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(histos,      cachefile, pickle.HIGHEST_PROTOCOL)


"""
loops over a set of files with common name to fill the migration matrices
"""        
def runUEAnalysis(inF,outF,wgtIdx,varIdx,cfgDir):

    print '[runAnalysis] %s -> %s wgtIdx=%d,varIdx=%d'%(inF,outF,wgtIdx,varIdx)
    
    #configure from pickle file
    ueHandler=UEAnalysisHandler(os.path.join(cfgDir,'analysiscfg.pck'))

    #loop over the tree to fill histos
    ue=UEEventCounter(EVENTAXES)
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
        gen_passSelPtTtbar = True if (t.gen_ptttbar>PTTTBAR_THR and gen_passSel) else False
        passSel=((t.passSel>>varIdx) & 0x1)
        passSelPtTtbar     = True if (t.ptttbar[varIdx]>PTTTBAR_THR and passSel) else False

        #fill histos (if sliceVar is None, non-sliced histos are filled)
        for sliceVar in SLICEQUANTILES.keys()+[None]:
            
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
            for obs in OBSQUANTILES:
                ueHandler.fillInclusive(sliceVarVals=sliceVarVals,
                                        obs=obs,
                                        ue=ue,
                                        weight=t.weight[wgtIdx],
                                        gen_passSel=gen_passSel,
                                        passSel=passSel)

                #loop over axes defining away/towards/transverse regions
                for a in EVENTAXES:
                    if a=='phittbar':
                        ueHandler.fillDifferential(sliceVarVals=sliceVarVals,
                                                   obs=obs,
                                                   a=a,
                                                   ue=ue,
                                                   weight=t.weight[wgtIdx],
                                                   gen_passSel=gen_passSelPtTtbar,
                                                   passSel=passSelPtTtbar)
                    else:
                        ueHandler.fillDifferential(sliceVarVals=sliceVarVals,
                                                   obs=obs,
                                                   a=a,
                                                   ue=ue,
                                                   weight=t.weight[wgtIdx],
                                                   gen_passSel=gen_passSel,
                                                   passSel=passSel)



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

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',    dest='input',  help='input',                       default='MC13TeV_TTJets_dilpowheg_0.root',   type='string')
    parser.add_option(      '--ptThr', dest='ptThr',  help='ptThreshold',                 default=0.9,   type=float)
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
        outDir=opt.out+'/analysis_%d_%d'%(opt.wgtIdx,opt.varIdx)
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
