import ROOT
import os
import numpy as np
from TopLJets2015.TopAnalysis.Plot import *
from PPSEfficiencyReader import PPSEfficiencyReader


baseDir='${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/'
xangle=120

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPaintTextFormat("4.2f");


for ppsEffFile in ['PreliminaryEfficiencies_October92019_1D2DMultiTrack.root', #'pixelEfficiencies.root'
               ]:
    ppsEffReader=PPSEfficiencyReader(os.path.join(baseDir,ppsEffFile))

    for rp in [23,123]:

        
        p=Plot('raddameff_%d'%rp,com='13 TeV')
        p.range=[0,1]
        p.xtit='#xi'
        p.ytit='Efficiency'
        for i,era in enumerate(['B','C','D','E','F']):
            gr=ROOT.TGraphErrors()
            gr.SetMarkerStyle(20+i)
            gr.SetTitle('2017'+era)
            for xi in np.linspace(0,0.2,50):
                ip=gr.GetN()
                eff,effUnc=ppsEffReader.getPPSEfficiency('2017'+era,xangle,xi,rp)
                gr.SetPoint(ip,xi,eff)
                gr.SetPointError(ip,0,effUnc)
            p.add(gr,'2017'+era, color=1, isData=False, spImpose=True, isSyst=False)
        p.show(outDir='./', lumi=37500)
