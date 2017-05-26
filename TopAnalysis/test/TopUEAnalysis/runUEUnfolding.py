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
    parser.add_option('-p', '--plotter', dest='plotter',  help='plotter with total distributions',        default=None,                    type='string')
    parser.add_option('-d', '--dir',     dest='dir',      help='chunk directory for pseudo-experiments',  default=None,                    type='string')
    parser.add_option('-t', '--tag',     dest='tag',      help='process this tag',                        default='MC13TeV_TTJets',        type='string')
    parser.add_option(      '--histo',   dest='histo',    help='histogram',                               default='chmult_None_inc',       type='string')
    parser.add_option(      '--sigName', dest='sigName',  help='signal name',                             default='t#bar{t}',              type='string')
    parser.add_option(      '--opt_tau', dest='opt_tau',  help='opt_tau',                                 default=-1,                      type=float)
    parser.add_option('-o', '--out',     dest='out',      help='output',                                  default='./UEanalysis/unfold',   type='string')
    (opt, args) = parser.parse_args()

    #prepare output
    os.system('mkdir -p %s'%opt.out)

    ROOT.gROOT.SetBatch(True)
  
    #pre-compile the unfolding script, if needed
    ROOT.gSystem.CompileMacro('${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/TopUEAnalysis/UEUnfold.C','k')
    from ROOT import UEUnfold

    #build the list of files
    file_list=[]
    if opt.dir.startswith('/store/') : 
        file_list = getEOSlslist(opt.dir)
    else:
        file_list=[os.path.join(opt.dir,x) for x in os.listdir(opt.dir) if re.match('%s_\d'%opt.tag,x)]

    #run unfolding
    opt_tau=opt.opt_tau
    fOut=ROOT.TFile.Open('%s/%s.root'%(opt.out,opt.histo),'RECREATE')
    globalBias,globalPulls=None,None
    binBias,binPulls=None,None
    for itoy in xrange(0,len(file_list)):
        f=file_list[itoy]
        results=UEUnfold(opt.histo,opt.plotter,f,opt.sigName, opt_tau,True if itoy==0 else False)

        fOutDir=fOut.mkdir('toy_%d'%itoy)
        for r in results: 
            rname=r.GetName()
            
            if rname in ['cscan','cunfolded','TVectorT<float>','corrected_data']:
                fOut.cd()
                if rname=='TVectorT<float>' :
                    if r.GetNoElements()==2: 
                        opt_tau=r[0]
                        rname='opt_tau'
                        print '[runUEUnfolding] Set optimal tau to ',opt_tau
                    else:
                        rname='chi2'
            else:
                fOutDir.cd()

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
                
            r.Write(rname)
        fOut.cd()
        results.Delete()

        #break

    print globalPulls.GetMean(),globalPulls.GetRMS()
    #all done
    fOut.cd()
    for x in [globalBias,globalPulls,binBias,binPulls]: x.Write()
    fOut.Close()



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
