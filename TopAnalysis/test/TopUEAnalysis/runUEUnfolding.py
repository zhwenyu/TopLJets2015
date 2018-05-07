#!/usr/bin/env/python

import sys
import os
import re
import optparse
import ROOT

def computeChisquare(data,model,cov):
    """ computes the chi^2 between two graphs and the covariance matrix associated to the first graph """

    #n-points
    np=data.GetN()
    
    #invert covariance matrix
    invCovM=None
    if cov:        
        nx,ny=cov.GetNbinsX(),cov.GetNbinsY()
        covM=ROOT.TMatrixF(nx,ny)
        for xbin in xrange(1,nx+1):
            for ybin in xrange(1,ny+1):
                covM[xbin-1][ybin-1]=cov.GetBinContent(xbin,ybin)
        invCovM=ROOT.TMatrixF(covM)
        invCovM.Invert()

    #vector of differences
    diffVec=ROOT.TVectorF(np)
    x,y=ROOT.Double(0),ROOT.Double(0)
    for i in xrange(0,np):

        data.GetPoint(i,x,y)
        ydata_i=float(y)

        model.GetPoint(i,x,y)
        ymodel_i=float(y)

        #the difference has to be multiplied by the bin width 
        ex=data.GetErrorX(i)
        diffVec[i]=ex*(ymodel_i-ydata_i)

    #chi^2 for distribution analysis
    chi2=0
    if invCovM:
        for i in xrange(0,np):
            for j in xrange(0,np):
                chi2 += diffVec[i]*invCovM[i][j]*diffVec[j]            
    else:
        for i in xrange(0,np):
            ey=data.GetErrorY(i)
            ex=data.GetErrorX(i)
            if ey*ex==0: continue
            chi2 += (diffVec[i]/(ey*ex))**2
        
    pval=ROOT.TMath.Prob(chi2,np-1)

    return chi2,np-1,pval


"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-s', '--step',        dest='step',         help='step [%default]',                                    default=1,                       type=int)
    parser.add_option(      '--plotter',     dest='plotter',      help='plotter with total distributions [%default]',        default=None,                    type='string')
    parser.add_option(      '--syst',        dest='systPlotter',  help='syst plotter [%default]',                            default=None,                    type='string')
    parser.add_option('-d', '--dir',         dest='dir',          help='chunk directory for pseudo-experiments [%default]',  default=None,                    type='string')
    parser.add_option('-t', '--tag',         dest='tag',          help='process this tag [%default]',                        default='MC13TeV_TTJets',        type='string')
    parser.add_option(      '--sigName',     dest='sigName',      help='signal name [%default]',                             default='t#bar{t}',              type='string')
    parser.add_option('-o', '--out',         dest='out',          help='output [%default]',                                  default='./UEanalysis/unfold',   type='string')
    (opt, args) = parser.parse_args()

    #prepare output
    os.system('mkdir -p %s'%opt.out)
    outName='%s/unfold_summary.root'%opt.out
  
    #pre-compile the unfolding script, if needed
    ROOT.gSystem.CompileMacro('${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/TopUEAnalysis/UEUnfold.C','k')
    unf=ROOT.UEUnfold()

    #unfold data
    if opt.step==1:
        print 'Unfolding data'
        unf.unfoldData(opt.plotter,opt.systPlotter,opt.sigName)
        results=unf.getResults()

        print 'Saving results in',outName
        fOut=ROOT.TFile.Open(outName,'RECREATE')
        for r in results: r.Write(r.GetName())
        fOut.Close()

    #run toys
    elif opt.step in [0,2]:

        #build the list of files
        file_list=[]
        if opt.dir.startswith('/store/') : 
            file_list = getEOSlslist(opt.dir)
        else:
            file_list=[os.path.join(opt.dir,x) for x in os.listdir(opt.dir) if re.match('%s_\d'%opt.tag,x)]

        chunkA = [ file_list[i] for i in xrange(0,len(file_list)+1,2) ]
        chunkB = [ file_list[i] for i in xrange(1,len(file_list),2) ]
        
        if opt.step==0:
            print 'Hadding 1/2 of the files for the migration matrix'
            os.system( 'hadd -f -k %s/ChunkAForToys.root %s'%(opt.out,' '.join( str(x) for x in chunkA) ) )

        if opt.step==2:

            #get a stat independent migration and fakes estimate
            fIn=ROOT.TFile.Open('%s/ChunkAForToys.root'%opt.out)
            indmig=fIn.Get('mig_0')
            indmig.SetDirectory(0)
            indfakes=fIn.Get('fakes_0')
            indfakes.SetDirectory(0)
            refRec=fIn.Get('reco_0')
            refRec.Add(indfakes,-1)
            refRec.SetDirectory(0)
            refGen=fIn.Get('gen')
            refGen.SetDirectory(0)
            fIn.Close()

            print 'Running %d toys'%len(chunkB)
            fOut=ROOT.TFile.Open(outName,'UPDATE')
            opt_tau=fOut.Get('TVectorT<float>;1')[0]
            datasub=fOut.Get('datasub')
            fakes=fOut.Get('fakes')
            mig=fOut.Get('migration')
            norm=fakes.Integral()+datasub.Integral()
            indfakes.Scale(fakes.Integral()/indfakes.Integral())
            indmig.Scale(mig.Integral()/indmig.Integral())

            pvalRatio=None
            globalBias,globalPulls=None,None
            binBias,binPulls=None,None
            for itoy in xrange(0,len(chunkB)):
                unf.reset();
                status=unf.unfoldToy(chunkB[itoy],opt_tau,indmig,norm,indfakes)
                if not status: continue            

                #bottom line test
                toySF=unf.curToyTruth_.Integral()/refGen.Integral()
                #refRec.Scale(unf.curToyRec_.Integral()/refRec.Integral())
                refRec.Scale(toySF)
                #refRec.Scale(unf.curFolded_unfolded_.Integral()/refRec.Integral())
                refGr=ROOT.TGraphErrors(refRec)
                toyRec=ROOT.TGraphErrors(unf.curToyRec_)                
                #toyRec=ROOT.TGraphErrors(unf.curFolded_unfolded_)
                smearedChi2=computeChisquare(toyRec,refGr,None)
                
                toyUnf=ROOT.TGraphErrors(unf.curToyUnf_)
                toyGen=ROOT.TGraphErrors(unf.curToyTruth_)
                refGen.Scale(toySF)
                refGenGr=ROOT.TGraphErrors(refGen)
                unfChi2=computeChisquare(toyUnf,refGenGr,unf.curCov_)

                toyRec.Delete()
                refGr.Delete()
                toyGen.Delete()
                toyUnf.Delete()
                refGenGr.Delete()

                if not pvalRatio:
                    pvalRatio=ROOT.TH1F('pvalRatio',';p-val(unf.)/p-val(smeared);',50,0,10)
                    pvalRatio.SetDirectory(0)
                if smearedChi2[2]>0: pvalRatio.Fill(unfChi2[2]/smearedChi2[2])

                fOut.cd()
                fOut.rmdir('toy_%d'%itoy)
                fOutDir=fOut.mkdir('toy_%d'%itoy)
                fOutDir.cd()

                for r in unf.getResults():
                
                    rname=r.GetName()

                    #bias summary
                    if rname=='toy_bias':
                        if globalBias is None:
                            maxVal=r.GetMaximum()*5
                            globalBias=ROOT.TH1F('global_bias',';Bias;',100,-maxVal,maxVal)
                            globalBias.SetDirectory(0)
                            binBias=ROOT.TH2F('bin_bias',';Bin number;Bias;',r.GetNbinsX(),0,r.GetNbinsX(),100,-maxVal,maxVal)
                            binBias.SetDirectory(0)

                        for xbin in xrange(1,r.GetNbinsX()+1):
                            xbias=r.GetBinContent(xbin)
                            globalBias.Fill(xbias)
                            binBias.Fill(xbin-1,xbias)

                    #pull summary
                    if rname=='toy_pull':
                        if globalPulls is None:                        
                            globalPulls=ROOT.TH1F('global_pulls',';Pull;',100,-5,5)
                            globalPulls.SetDirectory(0)
                            binPulls=ROOT.TH2F('bin_pulls',';Bin number;Pull;',r.GetNbinsX(),0,r.GetNbinsX(),100,-5,5)
                            binPulls.SetDirectory(0)
                        for xbin in xrange(1,r.GetNbinsX()+1):
                            xpull=r.GetBinContent(xbin)
                            globalPulls.Fill(xpull)
                            binPulls.Fill(xbin-1,xpull)
                    
                    r.Write()

            print 'global bias',globalBias.GetMean(),globalBias.GetRMS()
            print 'global pull',globalPulls.GetMean(),globalPulls.GetRMS()
            
            #all done
            fOut.cd()
            for x in [pvalRatio,globalBias,globalPulls,binBias,binPulls]: x.Write()
            fOut.Close()



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
