#!/usr/bin/env/python

import glob
import sys
import os
import optparse
import ROOT
import numpy as np
import array as array

from UETools import *

#GLOBAL VARIABLES TO DEFINE THE ANALYSIS

PTTTBAR_THR=150

VARTITLES={
    'ptttbar'  :'p_{T}(t#bar{t})',
    'phittbar' :'#phi(t#bar{t})',
    'ptpos'    :'p_{T}(l^+)',
    'phipos'   :'#phi(l^+)',
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
    varVals={
        'pt_ttbar' :[[],[]],
        'phi_ttbar':[[],[]],
        'mll'      :[[],[]],
        'dphill'   :[[],[]],
        'phill'    :[[],[]],
        'nj'       :[[],[]],
        'nch'      :[[],[]]
        }
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%100==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()
        if i>5000 : continue
        #reco level
        passSel=(t.passSel&0x1)
        gen_passSel=t.gen_passSel
        if not passSel or not gen_passSel: continue

        if t.gen_pt_ttbar>PTTTBAR_THR:
            varVals['pt_ttbar'][0].append( t.gen_pt_ttbar )
            varVals['pt_ttbar'][1].append( t.rec_pt_ttbar[0]-t.gen_pt_ttbar )
            varVals['phi_ttbar'][0].append(t.gen_phi_ttbar*180./ROOT.TMath.Pi())
            varVals['phi_ttbar'][1].append(ROOT.TVector2.Phi_mpi_pi(t.rec_phi_ttbar[0]-t.gen_phi_ttbar)*180./ROOT.TMath.Pi())

        varVals['mll'][0].append(t.gen_mll)
        varVals['mll'][1].append(t.mll-t.gen_mll)

        varVals['dphill'][0].append(t.gen_dphill*180./ROOT.TMath.Pi())
        varVals['dphill'][1].append(ROOT.TVector2.Phi_mpi_pi(t.dphill-t.gen_dphill)*180./ROOT.TMath.Pi())

        varVals['phill'][0].append(t.gen_phill*180./ROOT.TMath.Pi())
        varVals['phill'][1].append(ROOT.TVector2.Phi_mpi_pi(t.phill-t.gen_phill)*180./ROOT.TMath.Pi())

        varVals['nj'][0].append(t.gen_nj)
        varVals['nj'][1].append((t.nj[0]-2)-t.gen_nj)

        ue=UEEventCounter(t)
        varVals['nch'][0].append(sum(ue.gen_nch))
        varVals['nch'][1].append(sum(ue.rec_nch)-sum(ue.gen_nch))

    #prepare output
    fOut=ROOT.TFile.Open('%s/UEanalysis.root'%opt.out,'RECREATE')
    outDir=fOut.mkdir('sliceVars')

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
        h2d=ROOT.TH2F(var,';Gen. level;Resolution;%', len(genvarQ)-1,array.array('d',genvarQ),50,dvarQ[0],dvarQ[1])
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

        #save to ROOT file
        outDir.cd()
        h2d.SetDirectory(outDir)
        h2d.Write()
        resGr.Write()

    #compute correlations at generator level
    ROOT.gStyle.SetPaintTextFormat("4.0f")
    h2d=ROOT.TH2F('slicecorr',';Variable; Variable; Correlation (%)',len(varVals)-2,0,len(varVals)-2,len(varVals)-2,0,len(varVals)-2)
    xbin=0
    for var1 in varVals:
        if '_ttbar' in var1 : continue
        xbin+=1
        h2d.GetXaxis().SetBinLabel(xbin,VARTITLES[var1])
        ybin=0
        for var2 in varVals:
            if '_ttbar' in var2 : continue
            rho=np.corrcoef(varVals[var1][0],varVals[var2][0])[0][1]
            ybin+=1
            h2d.GetYaxis().SetBinLabel(ybin,VARTITLES[var2])
            h2d.SetBinContent(xbin,ybin,rho*100)
    c.Clear()
    h2d.Draw('colztext')
    tex=ROOT.TLatex()        
    tex.SetTextFont(42)
    tex.SetTextSize(0.035)
    tex.SetNDC()
    tex.DrawLatex(0.1,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.DrawLatex(0.8,0.96,'#sqrt{s}=13 TeV')
    c.Modified()
    c.Update()
    for ext in ['png','pdf']: c.SaveAs('%s/slicecorr.%s'%(opt.out,ext))

    #all done
    outDir.cd()
    h2d.SetDirectory(outDir)
    h2d.Write()
    fOut.Close()


"""
defines the variables to slice the UE measurement and saves their reco quantiles
"""
def defineAnalysisBinning(opt):


    fIn=ROOT.TFile.Open('%s/UEanalysis.root'%opt.out,'UPDATE')

    #
    # SLICE VARIABLES
    # readout resolution curves and determine the bins for the observables
    #
    print 'Defining axes for the slice variables'
    slicingAxes=[]
    for k in fIn.Get('sliceVars').GetListOfKeys():
        kname=k.GetName()
        if not '_resol' in kname: continue

        var=kname.split('_resol')[0]

        #special case for jet multiplicity
        if var=='nj':
            genBin=[0,1,2,3]
            recBin=[0,1,2,3]
        else:
             
            #get resolution map and quantiles       
            inigenBin=[]           
            genResol=[]
            resolGr=fIn.Get('sliceVars/%s'%kname)
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
        slicingAxes.append( (ROOT.TAxis(len(genBin)-1,array.array('d',genBin)),
                               ROOT.TAxis(len(recBin)-1,array.array('d',recBin)) ) )
        slicingAxes[-1][0].SetName('%s_genSlices'%var)
        slicingAxes[-1][1].SetName('%s_recSlices'%var)


    #
    # OBSERVABLES
    # use quantiles to determine the binning for the observables : nch, sumpt and avgpt
    #
    print 'Defining axes for the observables'
    obsVals={'nch'  :[[],[],[],[]], 
             'ptsum':[[],[],[],[]], 
             'avgpt':[[],[],[],[]]}
    t=ROOT.TChain('tue')
    for f in opt.input.split(','): t.AddFile(f)
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%100==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()
        if i>5000: break

        #gen level particle counting
        gen_passSel=t.gen_passSel
        if not gen_passSel: continue

        ue=UEEventCounter(t)
        
        for k in xrange(0,3):
            obsVals['nch'][k].append( ue.gen_nch[k] )
            obsVals['ptsum'][k].append( ue.gen_ptsum[k] )
            obsVals['avgpt'][k].append( ue.gen_avgpt[k] )
        
        incNch=sum(ue.gen_nch)
        incPtSum=sum(ue.gen_ptsum)
        incAvgPt=incPtSum/incNch if incNch>0 else 0.
        obsVals['nch'][3].append(incNch) 
        obsVals['ptsum'][3].append(incPtSum)
        obsVals['avgpt'][3].append(incAvgPt) 


    #determine quantiles and save as binnings
    obsAxes={}
    for obs in obsVals:
        obsAxes[obs]=[[],[]]
        for i in xrange(0,len(obsVals[obs])):
            regName=getRegionName(idx=i)
            genVarQ = np.percentile( np.array( obsVals[obs][i] ), OBSQUANTILES[obs] )
            genVarQ[0]=0.
            if i<3 : 
                obsAxes[obs][0].append( ROOT.TAxis(len(genVarQ)-1,array.array('d',genVarQ)) )
                obsAxes[obs][0][-1].SetName('%s%s_Obs'%(obs,regName))
            else :
                obsAxes[obs][1].append( ROOT.TAxis(len(genVarQ)-1,array.array('d',genVarQ)) )
                obsAxes[obs][1][-1].SetName('%s%s_Obs'%(obs,regName))



    #save binnings to ROOT file
    print 'Saving binnings'
    outDir=fIn.Get('bins')
    try:
        outDir.cd()
    except:
        outDir=fIn.mkdir('bins')
        outDir.cd()
    for axgen,axrec in slicingAxes:
        axgen.Write(axgen.GetName(),ROOT.TObject.kOverwrite)
        axrec.Write(axrec.GetName(),ROOT.TObject.kOverwrite)
    for obs in obsAxes:
        for i in xrange(0,len(obsAxes[obs])):
            for axobs in obsAxes[obs][i]:
                axobs.Write(axobs.GetName(),ROOT.TObject.kOverwrite)

    outDir=fIn.Get('analysisTemplates')
    try:
        outDir.cd()
    except:
        outDir=fIn.mkdir('analysisTemplates')
        outDir.cd()

    #define the migration matrices, gen level/rec level histos
    print 'Saving gen/rec level histos and migration matrix templates'
    #inclusive
    for obs in obsAxes:
        for i in xrange(0,len(obsAxes[obs])):

            name=obs
            if i==0: name += '_diff'
            nbins=0
            for axobs in obsAxes[obs][i]: nbins += axobs.GetNbins()
            for level in ['gen','rec']:           
                h=ROOT.TH1F('%s_%s'%(name,level),'%s_%s'%(name,level),nbins,0,nbins)
                h.SetDirectory(outDir)
                h.Write(h.GetName(),ROOT.TObject.kOverwrite)
            h=ROOT.TH2F('%s_m'%name,'%s_m'%name,nbins,0,nbins,nbins,0,nbins)
            h.SetDirectory(outDir)
            h.Write(h.GetName(),ROOT.TObject.kOverwrite)

    #sliced
    for axgen,axrec in slicingAxes:

        sliceName=axgen.GetName().split('_genSlices')[0]

        for obs in obsAxes:
            for i in xrange(0,len(obsAxes[obs])):

                name=sliceName+'_'+obs
                if i==0 : name += '_diff'
                nbinsObs=0
                for axobs in obsAxes[obs][i]: nbinsObs += axobs.GetNbins()
                nbinsX=axgen.GetNbins()*nbinsObs
                nbinsY=axrec.GetNbins()*nbinsObs

                for level,nbins in [('gen',nbinsX),('rec',nbinsY)]:           
                    h=ROOT.TH1F('%s_%s'%(name,level),'%s_%s'%(name,level),nbins,0,nbins)
                    h.SetDirectory(outDir)
                    h.Write(h.GetName(),ROOT.TObject.kOverwrite)
                h=ROOT.TH2F('%s_m'%name,'%s_m'%name,nbinsX,0,nbinsX,nbinsY,0,nbinsY)
                h.SetDirectory(outDir)
                h.Write(h.GetName(),ROOT.TObject.kOverwrite)


    #all done
    fIn.Close()


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
    if opt.step==2:

        #prepare output
        outDir=opt.out+'/analysis_%d_%d'%(opt.wgtIdx,opt.varIdx)
        os.system.mkdir(outDir)

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
            tasklist.append((filename,'%s/%s'%(outdir,baseFileName),opt.wgtIdx,opt.varIdx))

        if opt.queue=='local':
         if opt.jobs>1:
            print ' Submitting jobs in %d threads' % opt.jobs
            import multiprocessing as MP
            pool = MP.Pool(opt.jobs)
            pool.map(runUEAnalysisPacked,tasklist)
         else:
            for fileName,outfile,wgtIdx,varIdx in taskList:
                runUEAnalysis(fileName,outfile,wgtIdx,varIdx)
        else:
        cmsswBase=os.environ['CMSSW_BASE']
        for fileName,outfile,wgtIdx,varIdx in tasklist:
            localRun='python %s/src/TopLJets2015/TopAnalysis/test/TopUEAnalysis/runUEanalysis.py -i %s -o %s -q local -s 2 -w %d -v %d'%(cmsswBase,fileName,outfile,wgtIdx,varIdx)
            cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localRun)
            print cmd
            os.system(cmd)


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
