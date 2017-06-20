import ROOT

from os import listdir
from os.path import isfile, join
import re

fIn=ROOT.TFile.Open('UEanalysis/analysis/MC13TeV_TTJets.root')

systs=['','tkeff','tkeffeta','tkeffbcdef'] #,'tkeffgh']

colors={'':ROOT.kOrange-8,
        'tkeff':ROOT.kOrange-6,
        'tkeffbcdef':ROOT.kOrange-4,
        'tkeffgh':ROOT.kOrange+4,
        'tkeffeta':ROOT.kBlue+2}


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
from TopLJets2015.TopAnalysis.Plot import Plot

for dist in ['chmult_None_inc_syst_True',
             'chflux_None_inc_syst_True',
             'chfluxz_None_inc_syst_True',
             'chavgpt_None_inc_syst_True',
             'chavgpz_None_inc_syst_True',
             'aplanarity_None_inc_syst_True',
             'sphericity_None_inc_syst_True',
             'C_None_inc_syst_True',
             'D_None_inc_syst_True'
             ]:

    plots={}
    h=fIn.Get(dist)
    try :
        h.Integral()
    except:
        print 'Skipping ',dist
        continue

    for ybin in xrange(1,h.GetNbinsY()+1):
        label=h.GetYaxis().GetBinLabel(ybin)
        if not label in systs: continue
        plots[label]=h.ProjectionX('nominal' if label == '' else label,ybin,ybin)
        plots[label].SetDirectory(0)
        plots[label].SetTitle(label)

    p=Plot(dist)
    for h in plots:
        plots[h].Scale(1./plots[h].Integral())
        if h =='' : continue        
        p.add(plots[h],h,colors[h],False,False,False)
    p.add(plots[''],'nominal',1,True,True,False)
    p.ratiorange=(0.7,1.3)
    p.doPoissonErrorBars=False
    p.show(outDir=' ~/www/TopUE_ReReco2016/eras/',lumi=35900,noStack=True,noRatio=False)
    #raw_input()

fIn.Close()
