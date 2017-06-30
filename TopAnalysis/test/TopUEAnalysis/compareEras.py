import ROOT

from os import listdir
from os.path import isfile, join
import re

#list files
path='UEanalysis/analysis'
onlyfiles = [f for f in listdir(path) if 'MuonEG' in f]

dist='chmult_None_inc_None_True'
#dist='chflux_None_inc_None_True'
#dist='chfluxz_None_inc_None_True'
#dist='chavgpt_None_inc_None_True'
#dist='chavgpz_None_inc_None_True'
#dist='aplanarity_None_inc_None_True'
#dist='sphericity_None_inc_None_True'
#dist='C_None_inc_None_True'
#dist='D_None_inc_None_True'

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
p=Plot(dist)
colors={'B':ROOT.kOrange-8,
        'C':ROOT.kOrange-6,
        'D':ROOT.kOrange-4,
        'E':ROOT.kOrange+4,
        'F':ROOT.kBlue+2,
        'BCDEF':ROOT.kOrange-4,
        'G':ROOT.kBlue-4,
        'H':ROOT.kBlue-9,
        'GH':ROOT.kBlue-9,
        }
for era in plots:
    plots[era].Scale(1./plots[era].Integral())
    p.add(plots[era],era,colors[era],False,False,False)
total.Scale(1./total.Integral())
p.add(total,'2016',1,True,True,False)
p.ratiorange=(0.7,1.3)
p.doPoissonErrorBars=False
p.show(outDir=' ~/www/TopUE_ReReco2016/eras/',lumi=35900,noStack=True,noRatio=False) #True)
#raw_input()
