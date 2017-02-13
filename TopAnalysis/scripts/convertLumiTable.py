#!/usr/bin/env python

import os,sys
import json
import commands
import optparse
from TopLJets2015.TopAnalysis.miniAnalyzer_cfi import analysis
import ROOT

"""
runs lumicalc
"""
def getLumiTable(trig,opt,templateH=None):

    #run brilcalc
    extraOpt='--hltpath %s*'%trig if len(trig)>0 else ''
    rawLumiTable=commands.getstatusoutput('brilcalc lumi -b "STABLE BEAMS" --normtag %s -u /pb -i %s %s'%(opt.normtag,opt.lumiMask,extraOpt))[1].split('\n')

    #parse lumi output file
    valList=[]
    for line in rawLumiTable:
        tokens=line.split()
        try:
            run,lumisec=tokens[1].split(':')
            lumi=tokens[-2]
            valList.append( (int(run),float(lumi)) )
        except:
            pass

    #fill histogram
    nRuns=len(valList)
    try:
        h=templateH.Clone('lumisec%s'%trig)
        h.Reset('ICE')
        for i in xrange(0,nRuns):
            run,lumi=valList[i]
            binLabel='%d'%run
            for xbin in xrange(1,h.GetNbinsX()+1):
                ilabel=h.GetXaxis().GetBinLabel(xbin)
                if ilabel!=binLabel : continue
                h.SetBinContent(xbin,lumi)
                break
    except:
        h=ROOT.TH1F('lumisec%s'%trig,'lumisec;Run number;Lumi recorded (pb);',nRuns,0,nRuns)
        for i in xrange(0,nRuns):
            run,lumi=valList[i]
            h.GetXaxis().SetBinLabel(i+1,'%d'%run)
            h.SetBinContent(i+1,lumi)

    h.SetDirectory(0)
    return h


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-n', '--normtag',  dest='normtag',     help='normtag',       default='/afs/cern.ch/user/l/lumipro/public/Normtags/normtag_DATACERT.json',    type='string')
    parser.add_option('-o', '--out',      dest='out'   ,      help='output dir',    default='data/era2016/',    type='string')
    parser.add_option('-l', '--lumi',     dest='lumiMask',    help='json with list of good lumis', default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt')
    (opt, args) = parser.parse_args()

    #build the integ lumi tables
    lumiTables={}
    triggerPaths=analysis.triggersToUse
    lumiTables['']=getLumiTable(trig='',opt=opt)
    print 'Total lumi',lumiTables[''].Integral(),'/pb'
    for trig in triggerPaths:
        lumiTables[trig]=getLumiTable(trig=trig,opt=opt,templateH=lumiTables[''])
        print trig,lumiTables[trig].Integral(),'/pb'

    #open file to save lumitables
    fOut=ROOT.TFile.Open('%s/lumisec.root'%opt.out,'RECREATE')
    for key in lumiTables:
        print key,lumiTables[key].Integral(),'/pb'
        lumiTables[key].SetDirectory(fOut)
        lumiTables[key].Write('lumisec%s'%key)
    fOut.Close()

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
