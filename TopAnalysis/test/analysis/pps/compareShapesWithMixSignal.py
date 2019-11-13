#!/usr/bin/env python

import ROOT
import os
from TopLJets2015.TopAnalysis.Plot import *

mass='800'
histos=[]
for d in ['mixNone','mix'+mass]:

    print d

    t=ROOT.TChain('data')
    for f in os.listdir(d+'/Chunks/'):
        t.AddFile(os.path.join(d+'/Chunks/',f))

    t.Draw('mmiss >> h(40,0,2000)','wgt*(mixType==0 && mmiss>0 && cat==143)','goff')
    h=ROOT.gDirectory.Get('h')
    nevts=h.Integral()
    h.Reset('ICE')

    t.Draw('mmiss >> h(40,0,2000)','wgt*(mixType==1 && mmiss>0 && cat==143)','goff')
    h=ROOT.gDirectory.Get('h')
    h.Scale(nevts/h.Integral())
    title=d.replace('mix','')
    title=title.replace('None','data')
    histos.append( h.Clone('h'+title) )
    histos[-1].SetTitle(title)
    histos[-1].SetDirectory(0)
    h.Reset('ICE')


ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

COLOURS=[1,ROOT.kAzure+1,ROOT.kGreen+3,ROOT.kRed+1]

p=Plot('mix_mmiss_siginj'+mass,com='13 TeV, 2017B')
p.doPoissonErrorBars=False
p.ratiorange=(0.58,1.24)
p.spimposeWithErrors=True
for ic in [0,1]:
    title=histos[ic].GetTitle()
    histos[ic].GetXaxis().SetTitle('Missing mass [GeV]')
    histos[ic].GetYaxis().SetTitle('Events')

    if ic==0:
        p.add(h=histos[ic],
              title=title,
              color=COLOURS[ic],
              isData=True,
              spImpose=False,
              isSyst=False)
    else:
        for frac in [0.001,0.01,0.1]:
            h=histos[ic]
            idx=len(histos)-1
            histos.append(h.Clone('siginj%d_%d'%(ic,idx)))
            histos[-1].Scale(frac)
            histos[-1].Add(histos[0],1.-frac)
            title='f(sig)=%.g'%frac
            p.add(h=histos[-1],
                  title=title,
                  color=COLOURS[idx],
                  isData=False,
                  spImpose=False,
                  isSyst=False)
            
p.show(outDir='./',lumi=4776,noStack=True,saveTeX=False,noRatio=False)
p.reset()    
    
