#!/usr/bin/env/python

import sys
import os
import re
import optparse
import ROOT

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
    parser.add_option(      '--histo',       dest='histo',        help='histogram [%default]',                               default='chmult_None_inc',       type='string')
    parser.add_option(      '--sigName',     dest='sigName',      help='signal name [%default]',                             default='t#bar{t}',              type='string')
    parser.add_option('-o', '--out',         dest='out',          help='output [%default]',                                  default='./UEanalysis/unfold',   type='string')
    (opt, args) = parser.parse_args()

    #prepare output
    os.system('mkdir -p %s'%opt.out)
    outName='%s/%s.root'%(opt.out,opt.histo)
  
    #pre-compile the unfolding script, if needed
    ROOT.gSystem.CompileMacro('${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/TopUEAnalysis/UEUnfold.C','k')
    unf=ROOT.UEUnfold()

    #unfold data
    if opt.step==1:
        print 'Unfolding data'
        unf.unfoldData(opt.histo,opt.plotter,opt.systPlotter,opt.sigName)
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
            os.system( 'hadd -f %s/ChunkAForToys.root %s'%(opt.out,' '.join( str(x) for x in chunkA) ) )

        if opt.step==2:

            #get a stat independent migration and fakes estimate
            fIn=ROOT.TFile.Open('%s/ChunkAForToys.root'%opt.out)
            indmig=fIn.Get('%s_0_mig'%opt.histo)
            indmig.SetDirectory(0)
            indfakes=fIn.Get('%s_fakes_True'%opt.histo)
            indfakes.SetDirectory(0)
            fIn.Close()

            print 'Running %d toys'%len(chunkB)
            fOut=ROOT.TFile.Open('%s/%s.root'%(opt.out,opt.histo),'UPDATE')
            opt_tau=fOut.Get('TVectorT<float>;1')[0]
            datasub=fOut.Get('datasub')
            fakes=fOut.Get('fakes')
            mig=fOut.Get('migration')
            norm=fakes.Integral()+datasub.Integral()
            indfakes.Scale(fakes.Integral()/indfakes.Integral())
            indmig.Scale(mig.Integral()/indmig.Integral())

            globalBias,globalPulls=None,None
            binBias,binPulls=None,None
            for itoy in xrange(0,len(chunkB)):
                unf.reset();
                results=unf.unfoldToy(opt.histo,chunkB[itoy],opt_tau,indmig,norm,indfakes)
                
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
            for x in [globalBias,globalPulls,binBias,binPulls]: x.Write()
            fOut.Close()



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
