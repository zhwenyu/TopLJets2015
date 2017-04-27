#!/usr/bin/env python

import ROOT
import os
import sys

files={'MC':'root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/TagAndProbe2015/TnP_Z_pp5TeV_MC_Zmu10mu10_trk_v4.root',
       'Data':'root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/TagAndProbe2015/TnP_Z_pp5TeV_data_trk_v4.root'}

baseHistos={'mass':ROOT.TH1F('mass',';Invariant mass [GeV];Events',20,60,120),
        'pt':ROOT.TH1F('probept',';Probe transverse momentum [GeV];Events',20,0,100),
        'abseta':ROOT.TH1F('abseta',';Probe pseudo-rapidity;Events',20,0,2.1),
        'glbChi2':ROOT.TH1F('glbChi2',';Probe global #chi^{2};Events',20,0,15),
        'tkDxy':ROOT.TH1F('tkDxy',';Probe track dxy [cm];Events',20,-0.25,0.25),
        'tkDz':ROOT.TH1F('tkDz',';Probe track dz [cm]; Events',20,-0.6,0.6),
        'tkValidPixelHits':ROOT.TH1F('tkValidPixelHits',';Probe track valid pixel hits;Events',10,0,10),
        'combRelIsoPF04':ROOT.TH1F('combRelIsoPF04',';Probe PFRelIso(0.4);Events',20,0,0.2),
        'numberOfMatchedStations':ROOT.TH1F('numberOfMatchedStations',';Probe number of matched stations;Events',7,0,7),
        'glbValidMuHits':ROOT.TH1F('glbValidMuHits',';Probe number of valid muon hits;Events',50,0,50),
        'tkTrackerLay':ROOT.TH1F('tkTrackerLay',';Probe number of tracker layers;Events',20,0,20)
        }
for h in baseHistos:
    baseHistos[h].Sumw2()
    baseHistos[h].SetDirectory(0)

histos={}
for title in files:

    histos[title]={}
    for h in baseHistos:
        histos[title][h]=baseHistos[h].Clone('%s_%d'%(h,len(histos)))
        histos[title][h].SetDirectory(0)
        histos[title][h].SetTitle(title)
        if title=='Data':
            histos[title][h].SetFillStyle(0)
            histos[title][h].SetMarkerStyle(20)
        else:
            histos[title][h].SetFillStyle(1001)
            histos[title][h].SetFillColor(ROOT.kGray)

    fIn=ROOT.TFile.Open(files[title])
    t=fIn.Get('tpTree/fitter_tree')

    print 'Starting with %s with %d entries'%(title,t.GetEntriesFast())

    for i in xrange(0,t.GetEntriesFast()):
        
        t.GetEntry(i)
        
        #require the tag to have fired the trigger
        if t.tag_HIL2Mu15==0 : continue

        #require under the Z mass region
        if t.mass<60 or t.mass>120: continue

        #probe requirements
        probePassKin=True if t.abseta<2.1 and t.pt>18 else False
        probePassId=True  if t.TightHI==1 else False 
        probePassIso=True if t.combRelIsoPF04<0.15 else False
        
        #fill n-1
        for h in baseHistos:
            if h in ['abseta','pt'] and not (probePassId and probePassIso) : continue
            elif h in ['combRelIsoPF04'] and not (probePassKin and probePassId) : continue
            elif h in ['mass'] and not (probePassKin and probePassId and probePassIso) : continue
            elif not (probePassKin and probePassIso) : continue
            val=getattr(t,h)
            histos[title][h].Fill(val)
    fIn.Close()

import pickle
cachefile = open("muonhistos.pck", 'w')
pickle.dump(histos, cachefile, pickle.HIGHEST_PROTOCOL)
cachefile.close()

