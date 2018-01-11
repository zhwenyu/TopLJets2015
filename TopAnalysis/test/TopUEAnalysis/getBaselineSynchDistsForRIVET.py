import os
import pickle
import ROOT
import sys

plotList=[
    ('gen_ptll',       's01-x01-y01'),
    ('gen_ptlsum',     's01-x02-y01'),
    ('gen_nj',         's01-x03-y01'),
    ('gen_ptpos',      's01-x04-y01'),
    ('gen_mll',        's01-x05-y01'),
    ]

baseUrl='/afs/cern.ch/user/p/psilva/work/Top/CMSSW_8_0_28/src/TopLJets2015/TopAnalysis/UEanalysis_test/chmult/inc/ue_test.root'
yodaUrl='TOP-17-015_synch.yoda'

if len(sys.argv)>1: baseUrl=sys.argv[1]
if len(sys.argv)>2: yodaUrl=sys.argv[2]

fIn=ROOT.TFile.Open(baseUrl)
fOut=open(yodaUrl,'w')

for p,pname in plotList:

    #dump plot in yoda format
    h=fIn.Get(p)
    if not '_nj' in p: h.Rebin()
    h.Scale(1./h.Integral(0,h.GetXaxis().GetNbins()+1))

    fOut.write('BEGIN YODA_SCATTER2D /CMS_2017_TOP_17_015/%s\n'%pname)
    fOut.write('Path=/CMS_2017_TOP_17_015/%s\n'%pname)
    fOut.write('Type=Scatter2D\n')
    fOut.write('# xval xerr- xerr+ yval yerr- yerr+\n')
    for i in xrange(1,h.GetNbinsX()+1):
        x=h.GetXaxis().GetBinCenter(i)
        exlo=0.5*h.GetXaxis().GetBinWidth(i)
        exhi=exlo
        y=h.GetBinContent(i)
        eyhi=h.GetBinError(i)
        eylo=eyhi
        fOut.write('%.6g %.6g %.6g %.6g %.6g %.6g\n'%(x,exlo,exhi,y,eylo,eyhi))
    fOut.write('END YODA_SCATTER2D\n')
    fOut.write('\n')

fOut.write('\n\n')
fOut.close()
fIn.Close()
