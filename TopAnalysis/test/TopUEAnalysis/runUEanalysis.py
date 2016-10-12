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
from UETools import *

#GLOBAL VARIABLES TO DEFINE THE ANALYSIS

PTTTBAR_THR=150

VARTITLES={
    'ptttbar'  :'p_{T}(t#bar{t})',
    'phittbar' :'#phi(t#bar{t})',
    'ptpos'    :'p_{T}(l^{+})',
    'phipos'   :'#phi(l^{+})',
    'ptll'     :'p_{T}(l,l)',
    'phill'    :'#phi(ll)',
    'sumpt'    :'#Sigma p_{T}(l)',
    'mll'      :'M(l,l)',
    'dphill'   :'#Delta#phi(l,l)',
    'nj'       :'N(jets)',
    'chmult'   :'N(ch)',
    'chflux'   :'#Sigma p_{T}(ch)',
    'chavgpt'  :'#bar{p}_{T}(ch)'
    }

SLICEQUANTILES={  
    'ptttbar'  : [ 100*x/20 for x in xrange(0,21) ],
    'phittbar' : [ 100*x/20 for x in xrange(0,21) ],
    'ptpos'    : [ 100*x/20 for x in xrange(0,21) ],
    'phipos'   : [ 100*x/20 for x in xrange(0,21) ],
    'ptll'     : [ 100*x/20 for x in xrange(0,21) ],
    'phill'    : [ 100*x/20 for x in xrange(0,21) ],
    'sumpt'    : [ 100*x/20 for x in xrange(0,21) ],
    'mll'      : [ 100*x/20 for x in xrange(0,21) ],
    'dphill'   : [ 100*x/20 for x in xrange(0,21) ],
    'nj'       : [ 100*x/20 for x in xrange(0,21) ],
    'chmult'   : [ 100*x/20 for x in xrange(0,21) ]
    }

OBSQUANTILES={
    'chmult' :[100*x/10 for x in xrange(0,11)],
    'chflux' :[100*x/10 for x in xrange(0,11)],
    'chavgpt':[100*x/10 for x in xrange(0,11)]
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
    for var in SLICEQUANTILES : varVals[var]=[[],[]]
    ue=UEEventCounter()
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%100==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()

        if i>5000 : continue

        #require a pure event selected at reco and gen levels
        passSel=(t.passSel&0x1)
        gen_passSel=t.gen_passSel
        if not passSel or not gen_passSel: continue

        #fill resolution arrays
        for obs,isAngle in [ ('ptttbar',  False),
                             ('phittbar', True),
                             ('ptpos',    False),
                             ('phipos',   True),
                             ('ptll',     False),
                             ('phill',    True),
                             ('sumpt',    False),
                             ('mll',      False),
                             ('dphill',   True),
                             ('nj',       False)
                             ]:
            var=getattr(t,'gen_%s'%obs)
            dVar=getattr(t,obs)[0]-var
            if isAngle : 
                var  = var*180./ROOT.TMath.Pi()
                dVar = ROOT.TVector2.Phi_mpi_pi(dVar)*180./ROOT.TMath.Pi()
            varVals[obs][0].append( var )
            varVals[obs][1].append( dVar )

        #special case for charged particle counting
        ue.count(t)
        varVals['chmult'][0].append(ue.gen_chmult)
        varVals['chmult'][1].append(ue.rec_chmult-ue.gen_chmult)


    varResolutions={}

    #prepare canvas
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPalette(ROOT.kBird)
    c=ROOT.TCanvas('c','c',550,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.12)
    c.SetBottomMargin(0.1)

    #determine resolution for 10 quantiles
    for var in varVals:

        #use quantiles and extremes to define the 2D resolution histogram
        genvarQ  = np.percentile( np.array(varVals[var][0]), [100*x/10 for x in xrange(0,11)] )        
        if var=='nj': genvarQ = [0,1,2,5]

        dvarQ = np.percentile( np.array(varVals[var][1]), [2.5,97.5])
        h2d=ROOT.TH2F(var,
                      ';%s (gen. level);Resolution;'%VARTITLES[var]+'%', 
                      len(genvarQ)-1,array.array('d',genvarQ),50,dvarQ[0],dvarQ[1])
        for i in xrange(0,len(varVals[var][0])):
            h2d.Fill(varVals[var][0][i],varVals[var][1][i])
        
        #fit a gausian to determine resolution per quantile
        resGr=ROOT.TGraphErrors()
        resGr.SetName(var+'_resol')
        resGr.SetMarkerStyle(20)
        for xbin in xrange(1,h2d.GetNbinsX()+1):
            tmp=h2d.ProjectionY('tmp',xbin,xbin)
            tmp.Fit('gaus','MQ+')
            gaus=tmp.GetFunction('gaus')
            resGr.SetPoint(xbin-1,h2d.GetXaxis().GetBinCenter(xbin),gaus.GetParameter(1))
            resGr.SetPointError(xbin-1,0,gaus.GetParameter(2))

            #normalize entries in quantile
            total=tmp.Integral()
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
        tex.DrawLatex(0.1,0.96,'#bf{CMS} #it{simulation preliminary} %s'%VARTITLES[var])
        tex.DrawLatex(0.8,0.96,'#sqrt{s}=13 TeV')
        c.Modified()
        c.Update()
        for ext in ['png','pdf']: c.SaveAs('%s/%s_resol.%s'%(opt.out,var,ext))

        varResolutions[var]=(h2d,resGr)


    #compute correlations at generator level
    ROOT.gStyle.SetPaintTextFormat("4.0f")
    hcorr=ROOT.TH2F('slicecorr',';Variable; Variable; Correlation (%)',len(varVals),0,len(varVals),len(varVals),0,len(varVals))
    xbin=0
    for var1 in varVals:
        xbin+=1
        hcorr.GetXaxis().SetBinLabel(xbin,VARTITLES[var1])
        ybin=0
        for var2 in varVals:
            rho=np.corrcoef(varVals[var1][0],varVals[var2][0])[0][1]
            ybin+=1
            hcorr.GetYaxis().SetBinLabel(ybin,VARTITLES[var2])
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
            if var=='phill' or var=='dphill' or var=='phi_ttbar':             
                inigenBin=[i*18. for i in xrange(0,11)]
                maxRes=18.

            #define gen bin merging quantiles if needed according to the resolution
            i=0
            genBin=[ inigenBin[0] ]
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
    axes=['phittbar','phipos','phill']
    obsVals={}
    for obs in OBSQUANTILES: 
        obsVals[obs]={'inc':[[]]}
        for a in axes: obsVals[obs][a]=[[],[],[]]
    ue=UEEventCounter(axes)
    t=ROOT.TChain('tue')
    for f in opt.input.split(','): t.AddFile(f)
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%100==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()
        if i>500: break

        #require a pure event selected at reco and gen levels
        passSel=(t.passSel&0x1)
        gen_passSel=t.gen_passSel
        if not passSel or not gen_passSel: continue
        ue.count(t)

        obsVals['chmult']['inc'][0].append(ue.gen_chmult)
        obsVals['chflux']['inc'][0].append(ue.gen_chflux)
        obsVals['chavgpt']['inc'][0].append(ue.gen_chavgpt) 
        for a in axes:
            for k in xrange(0,3):
                obsVals['chmult'][a][k].append( ue.gen_chmult_wrtTo[a][k] )
                obsVals['chflux'][a][k].append( ue.gen_chflux_wrtTo[a][k] )
                obsVals['chavgpt'][a][k].append( ue.gen_chavgpt_wrtTo[a][k] )
        
    #determine quantiles and save as binnings
    obsAxes={}
    for obs in obsVals:
        for a in obsVals[obs]:
            for i in xrange(0,len(obsVals[obs][a])):
                genVarQ = np.percentile( np.array( obsVals[obs][a][i] ), OBSQUANTILES[obs] )
                genVarQ[0]=0.                
                obsAxes[ (obs,a,i) ] = ROOT.TAxis(len(genVarQ)-1,array.array('d',genVarQ))
                obsAxes[ (obs,a,i) ].SetName('%s%s_%d'%(obs,a,i))

    #
    # DEFINE MIGRATION MATRICES, GEN/REC LEVEL HISTOS
    #
    print 'Saving gen/rec level histos and migration matrix templates'
    histos={}    
    for obs in obsVals:
        for a in obsVals[obs]:

            nbins=0
            for i in xrange(0,len(obsVals[obs][a])): nbins+= obsAxes[ (obs,a,i) ].GetNbins()

            for level in [False,True]:
                levelStr='rec' if level else 'gen' 
                name='%s_%s%s'%(levelStr,a,obs)            
                histos[ (obs,a,level) ] = ROOT.TH1F(name,name,nbins,0,nbins)
                histos[ (obs,a,level) ].SetDirectory(0)

            name='m_%s%s'%(a,obs)
            histos[ (obs,a) ] = ROOT.TH2F(name,name,nbins,0,nbins,nbins,0,nbins)
            histos[ (obs,a) ].SetDirectory(0)

    #sliced
    for var,isRec in slicingAxes:
        if isRec: continue

        axX,         axY         = slicingAxes[(var,False)], slicingAxes[(var,True)]
        nslicebinsX, nslicebinsY = axX.GetNbins(),           axY.GetNbins()

        for obs in obsVals:
            for a in obsVals[obs]:

                nbinsObs=0
                for i in xrange(0,len(obsVals[obs][a])): nbinsObs += obsAxes[ (obs,a,i) ].GetNbins()

                nbinsX,nbinsY=nslicebinsX*nbinsObs, nslicebinsY*nbinsObs
                for level,nbins in [(False,nbinsX),(True,nbinsY)]:
                    levelStr='rec' if level else 'gen'
                    name='%s_%s%s_%s'%(levelStr,a,obs,var)
                    histos[ (obs,a,level,var) ] = ROOT.TH1F(name,name,nbins,0,nbins)
                    histos[ (obs,a,level,var) ].SetDirectory(0)

                name='m_%s%s_%s'%(a,obs,var)
                histos[ (obs,a,var) ] = ROOT.TH2F(name,name,nbinsX,0,nbinsX,nbinsY,0,nbinsY)
                histos[ (obs,a,var) ].SetDirectory(0)

    #all done, save to pickle file
    with open(os.path.join(opt.out,'analysiscfg.pck'), 'w') as cachefile:
        pickle.dump(obsAxes, cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(histos,  cachefile, pickle.HIGHEST_PROTOCOL)




"""
loops over a set of files with common name to fill the migration matrices
"""        
def runAnalysis(analysisCfg,inF,outF,wgtIdx,varIdx):

    #FIXME: implement read axis, histogram templates, filler function
    ueHandler=UEAnalysisHandler(analysisCfg)

    #loop over the tree
    t=ROOT.TChain('tue')
    t.AddFile(inF)
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%100==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()

        #selection flags
        gen_passSel        = t.gen_passSel
        gen_passSelPtTtbar = True if (t.gen_pt_ttbar>PTTTBAR_THR and gen_passSel) else False
        passSel=((t.passSel>>varIdx) & 0x1)

        #count particles
        ue=UEEventCounter(t)

        #FIXME: fill histograms using the ueHandler


    #save to ROOT file
    fOut=ROOT.TFile.Open(outF,'RECREATE')
    for h in ueHandler.histos:
        h.SetDirectory(fOut)
        h.Write()
    fOut.Close()

"""
Wrapper for when the analysis is run in parallel
"""
def runAnalysisPacked(args):
    try:
        fileName,outDir,varIdx,wgtIdx=args
        runUEAnalysis(fileName,outDir,varIdx,wgtIdx)
    except : # ReferenceError:
        print 50*'<'
        print "  Problem with", name, "continuing without"
        print 50*'<'
        return False


def main():

    eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',       dest='input',    help='input',                   default='MC13TeV_TTJets_dilpowheg_0.root',   type='string')
    parser.add_option('-s', '--step',     dest='step',     help='step',                    default=1,   type=int)
    parser.add_option('-w', '--wgt',      dest='wgtIdx',   help='weight index to use',     default=0,   type=int)
    parser.add_option('-v', '--var',      dest='varIdx',   help='calib index to use',      default=0,   type=int)
    parser.add_option('-j', '--jobs',     dest='jobs',     help='jobs to run in parallel', default=1,   type=int)
    parser.add_option('-o', '--out',      dest='out',      help='output',                  default='./UEanalysis',   type='string')
    (opt, args) = parser.parse_args()

    os.system('mkdir -p %s'%opt.out)

    if opt.step==0:
        determineSliceResolutions(opt)
    if opt.step==1:
        defineAnalysisBinning(opt)

#    if opt.step==2:
#
#        #prepare output
#        outDir=opt.out+'/analysis_%d_%d'%(opt.wgtIdx,opt.varIdx)
#        os.system.mkdir(outDir)
#
#        #create the tasklist
#        file_list=[]
#        if os.path.isdir(opt.input):
#            for file_path in os.listdir(opt.input):
#               if file_path.endswith('.root'):
#                    file_list.append(os.path.join(opt.input,file_path))
#        elif opt.input.startswith('/store/'):
#                    file_list = getEOSlslist(opt.input)
#        elif '.root' in opt.input:
#                    file_list.append(opt.input)
#
#        tasklist=[]
#        for filename in file_list:
#                baseFileName=os.path.basename(filename)
#                tag,ext=os.path.splitext(baseFileName)
#            if len(onlyList)>0:
#                    processThis=False
#                for filtTag in onlyList:
#                    if filtTag in tag:
#                       processThis=True
#            if not processThis : continue
#            tasklist.append((filename,'%s/%s'%(outdir,baseFileName),opt.wgtIdx,opt.varIdx))
#
#        if opt.queue=='local':
#         if opt.jobs>1:
#            print ' Submitting jobs in %d threads' % opt.jobs
#            import multiprocessing as MP
#            pool = MP.Pool(opt.jobs)
#            pool.map(runUEAnalysisPacked,tasklist)
#         else:
#            for fileName,outfile,wgtIdx,varIdx in taskList:
#                runUEAnalysis(fileName,outfile,wgtIdx,varIdx)
#        else:
#        cmsswBase=os.environ['CMSSW_BASE']
#        for fileName,outfile,wgtIdx,varIdx in tasklist:
#            localRun='python %s/src/TopLJets2015/TopAnalysis/test/TopUEAnalysis/runUEanalysis.py -i %s -o %s -q local -s 2 -w %d -v %d'%(cmsswBase,fileName,outfile,wgtIdx,varIdx)
#            cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localRun)
#            print cmd
#            os.system(cmd)
#

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
