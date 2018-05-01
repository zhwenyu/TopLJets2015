#!/usr/bin/env python

import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
#ROOT.gROOT.SetBatch(True)

c= ROOT.TCanvas("c","c",500,500)
c.SetTopMargin(0.06)
c.SetRightMargin(0.03)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.1)



fIn=ROOT.TFile.Open('/afs/cern.ch/work/p/psilva/TOP-17-010-final//analysis/Chunks/MC13TeV_TTJets_12.root')
for histo in ['mlb','ptlb','drlb','dphilb']:

    c.Clear()
    leg=ROOT.TLegend(0.15,0.93,0.6,0.69)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)

    histos=[]
    for cat,color,title in [('cor',ROOT.kRed,'correct'),('wro',1,'wrong')]:
        histos.append( fIn.Get(cat+histo).Clone() )
        histos[-1].SetTitle(title)
        histos[-1].SetLineWidth(2)
        histos[-1].SetLineColor(color)
        histos[-1].Scale(1./histos[-1].Integral(0,histos[-1].GetNbinsX()+1))
        leg.AddEntry(histos[-1],title,'l')
        histos[-1].GetYaxis().SetTitleOffset(1.4)
        histos[-1].GetYaxis().SetRangeUser(0,0.2)
        histos[-1].GetYaxis().SetTitle('PDF')
        if histo=='mlb': histos[-1].GetXaxis().SetRangeUser(20,200)
        histos[-1].Draw('hist' if len(histos)==1 else 'histsame')
    leg.Draw()

    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.DrawLatex(0.85,0.96,'13 TeV')

    c.Modified()
    c.Update()
    raw_input()
    for ext in ['png','pdf']:
        c.SaveAs('gen_%s.%s'%(histo,ext))
