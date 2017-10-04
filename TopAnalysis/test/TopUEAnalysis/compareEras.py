import ROOT

from os import listdir
from os.path import isfile, join
import re

#list files
path='store/TOP-17-015/chmult/inc'
onlyfiles = [f for f in listdir(path) if 'MuonEG' in f]

dist='reco_0'

plots={}
total=None
for f in onlyfiles:
    era=re.search('(?<=2016)\w+', f).group(0)[0]
    era='BCDEF' if era in ['B','C','D','E','F']  else 'GH'

    #getdistribution from file
    fIn=ROOT.TFile.Open(join(path,f))
    h=fIn.Get(dist)
    if not era in plots:
        plots[era]=h.Clone(era)
        plots[era].SetTitle(era)
        plots[era].SetDirectory(0)
    else:
        plots[era].Add(h)
    if total is None:
        total=h.Clone('total')
        total.SetTitle('2016')
        total.SetDirectory(0)
    else:
        total.Add(h)
    fIn.Close()

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
from TopLJets2015.TopAnalysis.Plot import Plot
p=Plot('chmult_eras')
colors={'B':ROOT.kRed+1,
        'C':ROOT.kRed+2,
        'D':ROOT.kRed+3,
        'E':ROOT.kRed+4,
        'F':ROOT.kRed+5,
        'BCDEF':ROOT.kRed,
        'G':ROOT.kCyan-1,
        'H':ROOT.kCyan-2,
        'GH':ROOT.kCyan,
        }
for era in plots:
    plots[era].Scale(1./plots[era].Integral())
    p.add(plots[era],era,colors[era],False,False,False)
total.Scale(1./total.Integral())
p.add(total,'2016',1,True,True,False)
p.ratiorange=(0.7,1.3)
p.cmsLabel='#bf{CMS} #it{preliminary}'
p.doPoissonErrorBars=False
p.ratioTitle='Ratio'
p.show(outDir=' ~/www/TOP-17-015/',lumi=35900,noStack=True,noRatio=False) #True)
#raw_input()
