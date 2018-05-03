#!/usr/bin/env python

import glob
import sys
import os
import optparse
import ROOT
import re
import numpy as np
import array as array
import pickle

from UEEventCounter import *
from UEAnalysisHandler import *
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist
from UETools import getPurStab,getNormalizedPerColumn,getConditionNumber,getBinForVariable


def printStandardLabel(extra=[]):
    """
    dummy function to print basic information
    """
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.1,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.DrawLatex(0.8,0.96,'#sqrt{s}=13 TeV')

    #print the extra
    if len(extra)==0 : return
    txt=ROOT.TPaveText(0.12,0.94,0.4,0.94-0.05*len(extra),'NDC')
    txt.SetTextFont(42)
    txt.SetTextAlign(12)
    txt.SetTextSize(0.035)
    txt.SetFillColor(0)
    txt.SetBorderSize(1)
    for e in extra: txt.AddText(e)
    return txt




def defineAnalysisBinning(opt,ptthreshold,cuts,outDir):
    """
    determines the resolutions needed to define the migration matrices
    """

    #build the chain
    t=ROOT.TChain('tue')
    for f in opt.input.split(','): 
        if len(f):
            t.AddFile(f)

    #loop the available events and readout the values of the variables to memory
    ue=UEEventCounter(ptthreshold=ptthreshold,cuts=cuts)
    totalEntries=t.GetEntries()
    varVals=[[],[],[]]
    fakeVals=[]
    failVals=[]
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%5000==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()

        #count particles in the event
        ue.count(t=t,debug=False,isMC=True)
        recVal=getattr(ue,'rec_%s'%opt.obs)[0] 
        val=getattr(ue,'gen_%s'%opt.obs) 

        #require a pure event selected at reco and gen levels
        if not ue.rec_passSel[0] : 
            if ue.gen_passSel: failVals.append(val)
            continue

        #fill fakes if gen selection has failed
        if not ue.gen_passSel and ue.rec_passSel[0]: 
            fakeVals.append(recVal)
            continue

        #fill resolution arrays
        deltaVal=recVal - val
        varVals[0].append( val )
        varVals[1].append( recVal )
        varVals[2].append( deltaVal )
        
#    ue.fout.cd()
#    ue.output_tuple.Write()
#    ue.fout.Close()

    #
    # GEN/RECO LEVEL BIN DEFINITION
    # INITIAL BINNING DEFINED BY QUANTILES (= STATS IN EVERY BIN)
    # FINAL BINNING IS SATURATED BY RESOLUTION
    #

    #use quantiles and extremes to define the 2D resolution histogram
    genvarQ = np.percentile( np.array(varVals[0]), [ 2*x for x in xrange(0,51) ] )
    recvarQ = np.percentile( np.array(varVals[1]), [ 2*x for x in xrange(0,51) ] )
    if opt.obs in ['chavgpt','chavgpz','chflux','chfluxz']:
        genvarQ=genvarQ[genvarQ>ptthreshold[0]]
        recvarQ=recvarQ[recvarQ>ptthreshold[1]]

    dvarQ   = np.percentile( np.array(varVals[2]), [2.5,97.5])
    resolH  = ROOT.TH2F('resol',';Generator level;#Delta(reco-gen);',len(genvarQ)-1,array.array('d',genvarQ),50,dvarQ[0],dvarQ[1])
    resolH.Sumw2()
    for i in xrange(0,len(varVals[0])): resolH.Fill(varVals[0][i],varVals[2][i])
        
    #fit a gausian to determine resolution per quantile
    resGr=ROOT.TGraphErrors()
    resGr.SetName('diffresol')
    resGr.SetMarkerStyle(20)
    profGr=resGr.Clone('diffprof')
    for xbin in xrange(1,resolH.GetNbinsX()+1):
            
        tmp=resolH.ProjectionY('tmp',xbin,xbin)
        total=tmp.Integral(0,tmp.GetNbinsX()+1)
        if total<10 : continue
                
        xcut=resolH.GetXaxis().GetBinLowEdge(xbin)
        mean,sigma=tmp.GetMean(),tmp.GetRMS()
        npts=profGr.GetN()
        profGr.SetPoint(npts,resolH.GetXaxis().GetBinCenter(xbin),mean)
        profGr.SetPointError(npts,0,sigma)
        resGr.SetPoint(npts,resolH.GetXaxis().GetBinCenter(xbin),sigma)
        resGr.SetPointError(npts,0,tmp.GetRMSError())

        tmp.Delete()

    #binning at gen level
    genBins=[min(genvarQ[0],0.)]
    if opt.obs in ['chavgpt','chavgpz']: genBins=[1.0]
    nsd=1.5
    lastAcceptResol=resGr.GetErrorY(0)
    for xgen_i in genvarQ:        
        if xgen_i<genBins[0] : continue
        dx=xgen_i-genBins[-1]
        if dx<0 : continue
        if dx<nsd*lastAcceptResol : continue
        genBins.append( xgen_i )
        lastAcceptResol=resGr.Eval(xgen_i)
    del genBins[1]

    #binning at reco level
    recBins=[]
    for i in xrange(0,len(genBins)-1):
        x_i=genBins[i]
        x_ii=0.5*(genBins[i+1]+genBins[i])
        if x_i<recvarQ[-1] and x_ii<recvarQ[-1]:
            recBins.append(x_i)
            recBins.append(x_ii)
        else:
            x_i=0.5*(recvarQ[-1]+recBins[-1])
            x_ii=recvarQ[-1]
            recBins.append(x_i)
            recBins.append(x_ii)
            break
    
    #final tweaks
    #if opt.obs in ['chavgpt','chavgpz']  : genBins[0],recBins[0]=0.9,0.9
    #if opt.obs in ['chflux','chfluxz']   : genBins[0],recBins[0]=0.9,0.9
    if opt.obs in ['C','D','sphericity'] : genBins[-1],recBins[-1]=1,1
    if opt.obs in ['aplanarity']         : genBins[-1],recBins[-1]=0.5,0.5

    #migration matrix and reco/gen level distributions
    migH  = ROOT.TH2F('mig',';Generator level;Reconstruction level;Events (a.u.)',len(genBins)-1,array.array('d',genBins),len(recBins)-1,array.array('d',recBins))
    migH.Sumw2()
    sqmigH  = ROOT.TH2F('sqmig',';Generator level;Reconstruction level;Events (a.u.)',len(genBins)-1,array.array('d',genBins),len(genBins)-1,array.array('d',genBins))
    sqmigH.Sumw2()
    recH = ROOT.TH1F('rec','reco;Variable; Events/bin (a.u.)', len(recBins)-1, array.array('d',recBins) )
    genH = ROOT.TH1F('gen','gen;Variable; Events/bin (a.u.)',  len(genBins)-1, array.array('d',genBins) )
    failH = ROOT.TH1F('fail','fail;Variable; Events/bin (a.u.)',  len(genBins)-1, array.array('d',genBins) )
    for i in xrange(0,len(varVals[0])):
        genVal,recVal=varVals[0][i],varVals[1][i]

        migH.Fill(genVal,recVal)
        sqmigH.Fill(genVal,recVal)

        recBin=getBinForVariable(recVal,recH)
        recH.Fill(recVal,1./recH.GetXaxis().GetBinWidth(recBin))

        genBin=getBinForVariable(genVal,genH)
        genH.Fill(genVal,1./genH.GetXaxis().GetBinWidth(genBin))
        
    #distribution for failing reco step
    for x in failVals:
        failBin=getBinForVariable(x,failH)
        failH.Fill(x,1./failH.GetXaxis().GetBinWidth(failBin))
    effGr=ROOT.TGraph()
    effGr.SetLineWidth(2)
    effGr.SetLineColor(ROOT.kCyan-3)
    for xbin in xrange(1,genH.GetNbinsX()+1):
        totalFail=failH.GetBinContent(xbin)
        total=genH.GetBinContent(xbin)+totalFail
        if total==0 : continue
        n=effGr.GetN()
        effGr.SetPoint(n,n+1,totalFail/total)

    #contribution from fakes
    fakesH = ROOT.TH1F('fakes','fakes;Variable; Events/bin (a.u.)', len(recBins)-1, array.array('d',recBins) )
    for x in fakeVals: 
        recBin=getBinForVariable(x,fakesH)
        fakesH.Fill(x,1./fakesH.GetXaxis().GetBinWidth(recBin))
    fakesGr=ROOT.TGraph()
    fakesGr.SetLineWidth(2)
    fakesGr.SetLineStyle(9)
    fakesGr.SetLineColor(ROOT.kGray)
    for xbin in xrange(1,fakesH.GetNbinsX()+1,2):
        totalFakes=fakesH.Integral(xbin,xbin+1,"width")
        total=recH.Integral(xbin,xbin+1,"width")+totalFakes
        if total==0 : continue
        n=fakesGr.GetN()
        fakesGr.SetPoint(n,n+1,totalFakes/total)

    #normalized,inclusive versions of the 2D plots + purity/stability graphs
    normResolH=getNormalizedPerColumn(resolH)
    incResolH=resolH.ProjectionY('incResolH',1,resolH.GetNbinsX())
    incResolH.Scale(1./incResolH.Integral())
    incResolH.SetMarkerStyle(20)
    incResolH.GetYaxis().SetTitle('PDF')
    normMigH=getNormalizedPerColumn(migH)    
    normSqMigH=getNormalizedPerColumn(sqmigH)    
    purGr,stabGr,_=getPurStab(sqmigH)
    condK=getConditionNumber(migH)

    cutText=[]
    for cutKey in cuts:
        if cutKey=='region': cutText.append('%s region=%s'%(VARTITLES[cuts[cutKey][0]],cuts[cutKey][1]) )
        else               : cutText.append('%3.1f#leq%s<%3.1f'%(cuts[cutKey][0],VARTITLES[cutKey],cuts[cutKey][1]))

    condKText=[]
    if condK:
        condKText.append('#lambda_{max}=%3.1f'%condK[0])
        condKText.append('#lambda_{min}=%3.1f'%condK[1])
        condKText.append('cond K=%.6g'%condK[2])


    #show plot
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
    ROOT.gStyle.SetPaintTextFormat("3.0f");
    c=ROOT.TCanvas('c','c',550,500)
    c.SetTopMargin(0.06)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    c.Clear()
    allTexts=[]

    c.SetRightMargin(0.12)
    for h in [resolH,normResolH,migH,normMigH,sqmigH,normSqMigH]:

        drawOpt='colz'
        if h in [normResolH,normMigH,normSqMigH] : drawOpt+='text'
        h.Draw(drawOpt)
        h.GetYaxis().SetTitleOffset(1.1)

        extraText=['#bf{%s}'%VARTITLES[opt.obs]]+cutText[:]
        if h in [resolH,normResolH]: 
            profGr.Draw('p')
        if h in [migH,normMigH,sqmigH,normSqMigH]: 
            extraText+=condKText
        allTexts.append( printStandardLabel(extraText) )
        allTexts[-1].Draw()
        
        c.RedrawAxis()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']: c.SaveAs('%s/%s.%s'%(outDir,h.GetName(),ext))

    #inclusive resolution
    c.Clear()
    c.SetRightMargin(0.02)

    for p in [incResolH,
              resGr,
              (recH,genH,fakesH),
              (purGr,stabGr,effGr,fakesGr)]:
        cname=''
        extraText=['#bf{%s}'%VARTITLES[opt.obs]]+cutText[:]
        leg=ROOT.TLegend(0.65,0.94,0.98,0.8)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        if p==resGr:
            cname='resdiff'
            resGr.Draw('aep')
            resGr.GetXaxis().SetTitle('Generator level')
            resGr.GetYaxis().SetTitle('#Delta(reco-gen)')
            leg.AddEntry(p,'diff. resolution','ep')
        if p==incResolH:
            cname='incresol'
            incResolH.Draw("ep")
            incResolH.GetYaxis().SetTitleOffset(1.2)
            mean,sigma=incResolH.GetMean(),incResolH.GetRMS()
            incResolH.Fit('gaus','MRQ+','same',mean-sigma,mean+sigma)
            gaus=incResolH.GetFunction('gaus')
            mean=gaus.GetParameter(1)
            sigma=gaus.GetParameter(2)
            leg.AddEntry(incResolH,'inc. resolution','ep')
            leg.AddEntry(gaus,'gaussian fit','l')
            extraText.append('#mu=%3.2f #sigma=%3.2f'%(mean,sigma))
        if p==(recH,genH,fakesH):
            cname='dist'
            p[0].Draw('hist')
            p[0].GetYaxis().SetTitleOffset(1.1)
            p[0].SetLineWidth(2)
            p[1].Draw('histsame')
            p[1].SetLineWidth(2)
            p[1].SetLineColor(2)
            p[2].Draw('histsame')
            p[2].SetLineWidth(2)
            p[2].SetLineColor(ROOT.kGray)
            p[2].SetLineStyle(2)
            leg.AddEntry(p[0],'Reco','l')
            leg.AddEntry(p[1],'Gen','l')
            leg.AddEntry(p[2],'Fakes','l')
        if p==(purGr,stabGr,effGr,fakesGr):
            cname='purstab'
            purGr.Draw('al')
            purGr.GetXaxis().SetTitle('Bin number')
            purGr.GetYaxis().SetTitle('Purity, Stability, Efficiency or Fake fraction')
            purGr.GetYaxis().SetRangeUser(0,1)
            stabGr.Draw('l')
            effGr.Draw('l')
            fakesGr.Draw('l')
            leg.AddEntry(purGr,'purity','l')
            leg.AddEntry(stabGr,'stability','l')
            leg.AddEntry(effGr,'efficiency','l')
            leg.AddEntry(fakesGr,'fakes','l')
        
        leg.Draw()
        allTexts.append( printStandardLabel(extraText))
        allTexts[-1].Draw()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']: c.SaveAs('%s/%s.%s'%(outDir,cname,ext))


    #all done, save to pickle file
    analysisCfg={
        ('gen','axis')   : ROOT.TAxis(len(genBins)-1,array.array('d',genBins)),
        ('gen','histo')  : ROOT.TH1F('gencount',';Bin number;Events',len(genBins)-1,0,len(genBins)-1),
        ('reco','axis')  : ROOT.TAxis(len(recBins)-1,array.array('d',recBins)),
        ('reco','histo') : ROOT.TH1F('recocount',';Bin number;Events',len(recBins)-1,0,len(recBins)-1),
        ('mig','histo')  : ROOT.TH2F('migcount',';Generator level bin;Reconstruction level bin;Events',len(genBins)-1,0,len(genBins)-1,len(recBins)-1,0,len(recBins)-1)
        }
    plotCollection={
        'resol':resolH,
        'resol_norm':normResolH,
        'mig':migH,
        'mig_norm':normMigH,
        'sqmig':sqmigH,
        'sqmig_norm':normSqMigH,
        'resol_inc':incResolH,
        'resol_prof':profGr,
        'resol_diff':resGr,
        'rec':recH,
        'gen':genH,
        'purity':purGr,
        'stability':stabGr
        }
    
    with open(os.path.join(outDir,'analysiscfg.pck'), 'w') as cachefile:
        pickle.dump(analysisCfg,    cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(cuts,           cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(opt.obs,        cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(ptthreshold,    cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(plotCollection, cachefile, pickle.HIGHEST_PROTOCOL)

"""
loops over a set of files with common name to fill the migration matrices
"""        
def runUEAnalysis(inF,outF,cfgDir):

    print '[runAnalysis] %s -> %s (cfg @ %s)'%(inF,outF,cfgDir)
    
    isData=True if 'Data13TeV' in inF else False
    isTTJets=True if 'MC13TeV_TTJets' in inF else False
    ueHandler=UEAnalysisHandler(os.path.join(cfgDir,'analysiscfg.pck'),isTTJets)

    #loop over the tree to fill histos
    systList=SYSTS if isTTJets else [('',   0,0,False)]
    ue=UEEventCounter(ptthreshold=ueHandler.ptthreshold,
                      cuts=ueHandler.cuts,
                      systList=systList
                      )
    t=ROOT.TChain('tue')
    t.AddFile(inF)
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):

        t.GetEntry(i)

        if i%5000==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()

        #count particles
        ue.count(t=t,isMC=isTTJets)
        if isData: ue.w[0]=1.0
        ueHandler.fillHistos(ue)

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
        fileName,outfile,cfgDir=args
        runUEAnalysis(fileName,outfile,cfgDir)
    except : 
        print 50*'<'
        print "  Problem with", name, "continuing without"
        print 50*'<'
        return False


def main():

    ROOT.gROOT.SetBatch(True)

    #readout the configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',    dest='input',  help='input',                       default='MC13TeV_TTJets_dilpowheg_0.root',   type='string')
    parser.add_option(      '--ptThr', dest='ptThr',  help='ptThreshold gen,reco',        default='1.0,0.9',   type='string')
    parser.add_option('-s', '--step',  dest='step',   help='step',                        default=1,   type=int)
    parser.add_option(      '--obs',   dest='obs',    help='observable',                  default='chmult', type='string')
    parser.add_option(      '--slice', dest='slice',  help='slice',                       default=None,     type='string')
    parser.add_option(      '--reg',   dest='reg',    help='region',                      default=None,     type='string')
    parser.add_option('-j', '--jobs',  dest='jobs',   help='jobs to run in parallel',     default=1,   type=int)
    parser.add_option(      '--dryRun',  dest='dryRun',   help='do not submit',     default=False,  action='store_true')
    parser.add_option('-o', '--out',   dest='out',    help='output',                      default='./UEanalysis',   type='string')
    parser.add_option(      '--only',  dest='only',   help='csv list of tags to process', default='',  type='string')
    parser.add_option(      '--onlyExact',  dest='onlyExact',   help='set to true if only tag is to be fully matched', default=False,  action='store_true')
    parser.add_option('-q', '--queue', dest='queue',  help='Batch queue to use [default: %default]', default='local')
    (opt, args) = parser.parse_args()


    #define how the analysis should be configured
    if opt.step==1: 
        
        #prepare output and cuts
        cuts={}
        outDir='%s/%s'%(opt.out,opt.obs)
        if opt.slice : 
            outDir += '/%s'%opt.slice

            var="".join(re.findall("[a-zA-Z]",opt.slice))
            cuts[var]=[float(x) for x in re.findall("[-+]?\d*\.*\d+",opt.slice)]
        else:
            outDir += '/inc'
        if opt.reg :
            outDir += '_%s'%opt.reg
            regAxis,reg=opt.reg.split('=')
            regIdx=0 
            if reg=='tra':regIdx=1
            if reg=='awa':regIdx=2
            cuts['region']=(regAxis,regIdx)
        os.system('mkdir -p %s'%outDir)
        print 'Will apply the following cuts',cuts

        ptthreshold=[float(x) for x in opt.ptThr.split(',')]
    
        defineAnalysisBinning(opt,ptthreshold,cuts,outDir)
        
        with open('lastUE.dat','w') as logF:
            logF.write(outDir)


    #run the analysis
    elif opt.step==2:

        #tags to process
        onlyList=opt.only.split(',')

        #prepare output
        os.system('mkdir -p %s/Chunks'%opt.out)

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
                    if opt.onlyExact:
                        if tag.replace(filtTag+'_','').isdigit():
                            processThis=True
                    else:
                        if filtTag in tag:
                            processThis=True
                if not processThis : continue
            tasklist.append((filename,'%s/Chunks/%s'%(opt.out,baseFileName),opt.out))

        #run jobs locally
        if opt.queue=='local':
            if opt.jobs>1:
                print ' Running %d jobs in %d threads' % (len(tasklist), opt.jobs)
                import multiprocessing as MP
                pool = MP.Pool(opt.jobs)
                pool.map(runUEAnalysisPacked,tasklist)
            else:
                for fileName,outfile,cfgDir in tasklist:
                    runUEAnalysis(fileName,outfile,cfgDir)

        #submit jobs
        else:
            print 'Running %d jobs to the batch'%len(tasklist)
            cmsswBase=os.environ['CMSSW_BASE']
            FarmDirectory='{0}/src/TopLJets2015/TopAnalysis/{1}/FARM-UEANA'.format(cmsswBase,opt.out)
            os.system('mkdir -p %s'%FarmDirectory)
            print 'Scripts and logs will be available @ %s'%FarmDirectory

            condorScript='%s/condor.sub'%opt.out
            with open (condorScript,'w') as condor:

                condor.write('executable = {0}/$(jobName).sh\n'.format(FarmDirectory))
                #condor.write('output     = {0}/output_$(jobName).out\n'.format(FarmDirectory))
                #condor.write('error      = {0}/output_$(jobName).err\n'.format(FarmDirectory))
                condor.write('output = /dev/null\n')
                condor.write('error  = /dev/null\n')
                condor.write('log    = /dev/null\n')
                condor.write('+JobFlavour = "{0}"\n'.format(opt.queue))

                for fileName,_,_ in tasklist:
                    
                    jobName='%s'%(os.path.splitext(os.path.basename(fileName))[0])

                    jobScript='%s/%s.sh'%(FarmDirectory,jobName)
                    with open(jobScript,'w') as job:
                        job.write('#!/bin/bash\n')
                        job.write('WORKDIR=`pwd`\n')
                        job.write('echo "Working directory is ${WORKDIR}"\n')
                        job.write('cd %s\n'%cmsswBase)
                        job.write('eval `scram r -sh`\n')
                        job.write('cd ${WORKDIR}\n')
                        job.write('python {0}/src/TopLJets2015/TopAnalysis/test/TopUEAnalysis/runUEanalysis.py -q local -s 2 -i {1} -o {0}/src/TopLJets2015/TopAnalysis/{2}\n'.format(cmsswBase,fileName,opt.out))
                        job.write('echo "All done"\n')

                    os.system('chmod u+x %s'%jobScript)
                    condor.write('jobName=%s\n'%jobName)
                    condor.write('queue 1\n')

            if not opt.dryRun:
                os.system('condor_submit %s'%condorScript)


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
