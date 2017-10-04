import ROOT

from os import listdir
from os.path import isfile, join
import re

fIn=ROOT.TFile.Open('store/TOP-17-015/chmult/inc/MC13TeV_TTJets.root')

systs=[('nominal','_0',1,True,False),
       #('SF_{vtx.} (#mu) ','_22',ROOT.kGray,False,False),
       ('SF_{#eta} (from #mu)', '_25',ROOT.kRed,False,False),
       ('SF_{#eta} (from D*)','_26',ROOT.kCyan,False,False)] 


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
from TopLJets2015.TopAnalysis.Plot import Plot


p=Plot('chmult_tkeffsyst')
p.ratiorange=(0.7,1.3)
p.doPoissonErrorBars=False
p.cmsLabel='#bf{CMS} #it{simulation}'
p.doMCOverData = False
p.ratioTitle='Ratio '
for s,pf,color,isref,spimpose in systs:
    h=fIn.Get('reco'+pf)
    try :
        h.Scale(1./h.Integral())
    except:
        print 'Skipping ',s
        continue

    p.add(h,s,color,isref,spimpose,False)
    

p.show(outDir=' ~/www/TOP-17-015',lumi=35900,noStack=True,noRatio=False)

fIn.Close()
