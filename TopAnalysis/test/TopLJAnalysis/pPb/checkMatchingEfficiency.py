#!/usr/bin/env python

import ROOT
import sys

fIn=ROOT.TFile.Open(sys.argv[1])
passedBins=[5,3,9,10,11]
total=ROOT.TH1F('total','total',5,0,5)
passed=ROOT.TH1F('passed','passed',5,0,5)

#start table
print '\\hline'
print '\\multirow{2}{*}{Category} & \\multicolumn{2}{c}{Kinematics} & \\multicolumn{3}{c}{Algorithm} \\\\'
print '                          & 2b & 2q                        & $W_{\\rm had}$ & $\cPqt_{\\rm had}$ & $\cPqt_{\\rm lep}$ \\\\'
print '\\hline'

#fill table
for cat in ['1l4j2q','1l4j1b1q','1l4j2b']:
    print cat,fIn
    h=fIn.Get('jeteff_%s'%cat)
    for i in xrange(1,6):
        total.SetBinContent(i,h.GetBinContent(1))
        total.SetBinError(i,h.GetBinError(1))
        passed.SetBinContent(i,h.GetBinContent( passedBins[i-1] ) )
        passed.SetBinError(i,h.GetBinError( passedBins[i-1] ) )
    eff=ROOT.TEfficiency(passed,total)
    print '%10s'%cat,
    for k in xrange(1,6):
        print '& $%3.3f\\pm%3.3f$'%(eff.GetEfficiency(k),
                                    0.5*(eff.GetEfficiencyErrorUp(k)+eff.GetEfficiencyErrorLow(k))),
    print '\\\\'

    eff.Delete()

#all done
print '\\hline'
