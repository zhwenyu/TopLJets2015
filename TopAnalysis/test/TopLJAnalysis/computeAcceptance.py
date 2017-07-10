#!/usr/bin/env python

import ROOT
import sys
import numpy as np

fIn=ROOT.TFile.Open(sys.argv[1])

norm=fIn.Get('wgtcounter')
fid=fIn.Get('fidcounter')
fid.Divide(norm)

acc=None
uncs={
    'scale':[],
    'pdf':[],
    'aS':[],
    'npdf':[]
    }
for i in xrange(0,fid.GetNbinsX()):
    xbin=i+1
    val,unc=fid.GetBinContent(xbin),fid.GetBinError(xbin)
    if i==0: acc=(val,unc)
    if i in [1,2,3,4,6,8]: 
        uncs['scale'].append(val-acc[0])
    if i>=9 and i<=108:
        uncs['pdf'].append(val-acc[0])
    if i==109 or i==110:
        uncs['aS'].append(val-acc[0])
    if i>=284 and i<=324:
        uncs['npdf'].append(val-acc[0])

print 'Acceptance estimate'
scale=max(abs(max(uncs['scale'])),abs(min(uncs['scale'])))
pdf=np.sqrt(np.mean(np.array(uncs['pdf'])**2))
pdf=np.sqrt(pdf**2+(0.5*(max(uncs['aS'])+min(uncs['aS'])))**2)
npdf=np.sqrt(np.mean(np.array(uncs['npdf'])**2))
total=np.sqrt(acc[1]**2+scale**2+pdf**2+npdf**2)

print '%3.4f +/- %3.4f (stat) +/- %3.4f (scales) +/- %3.4f (PDF) +/- %3.4f (nPDF)'%(acc[0], acc[1], scale,pdf,npdf)
print '%3.4f +/- %3.4f (total)'%(acc[0],total)




