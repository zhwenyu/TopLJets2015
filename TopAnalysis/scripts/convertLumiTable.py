#!/usr/bin/env python

import ROOT
import os
import sys

inFile=sys.argv[1]

#parse lumi output file
valList=[]
for line in open(sys.argv[1], 'r'):
     tokens=line.split()
     try:
          run,lumisec=tokens[1].split(':')
          lumi=tokens[-2]
          valList.append( (int(run),float(lumi)/1.0e+06) )
     except:
          pass

#fill histogram
nRuns=len(valList)
h=ROOT.TH1F('lumisec','lumisec;Run number;Lumi recorded (pb);',nRuns,0,nRuns)
for i in xrange(0,nRuns):
     run,lumi=valList[i]
     h.GetXaxis().SetBinLabel(i+1,'%d'%run)
     h.SetBinContent(i+1,lumi)
totalLumi=h.Integral()
outFile='%s/lumisec.root'%os.path.dirname(inFile)
h.SaveAs(outFile)

print 'Histogram with lumi per run can be found in',outFile
print 'Total lumi=%3.1f pb'%totalLumi
