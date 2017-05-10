#!/usr/bin/env python

import ROOT
import pickle

cachefile = open("muonhistos.pck", 'r')
histos = pickle.load(cachefile)
cachefile.close()

c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetRightMargin(0.05)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.1)
c.SetLogy()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
for h in histos['Data']:
    mcScale=histos['Data'][h].Integral(5,15)/histos['MC'][h].Integral(5,15)

    histos['MC'][h].Scale(mcScale)
    histos['MC'][h].Draw('hist')
    histos['Data'][h].Draw('same')
    histos['MC'][h].GetYaxis().SetTitleOffset(1.2)
    #histos['MC'][h].GetYaxis().SetRangeUser(0.5,1.2*histos['MC'][h].GetMaximum())

    leg=c.BuildLegend(0.75,0.87,0.9,0.75)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    #leg.SetNColumns(2)
    leg.SetTextSize(0.035)
    txt=ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.035)
    txt.SetNDC()
    txt.DrawLatex(0.15,0.9,'#bf{CMS} #it{preliminary}')
    txt.DrawLatex(0.7,0.9,'25.8 pb^{-1} (5 TeV)')
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s.%s'%(h,ext))
        
