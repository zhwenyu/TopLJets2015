#!/usr/bin/env python

import ROOT

PLOTS=[]
for p in ['pt_l','y_l','pt_b','y_b','pt_wl','y_wl','pt_wjj','y_wjj']:
    for c in ['1l4j2b','1l4j1b1q','1l4j2q']:
        PLOTS.append(p+'_'+c)

mc=[ ('plots/MC8.16TeV_TTbar_pPb_mu/controlplots.root',353),
     ('plots/MC8.16TeV_TTbar_pPb_e/controlplots.root',292)]
data=[ 'plots/Data8.16TeV_pPb_mu/controlplots.root','plots/Data8.16TeV_pPb_e/controlplots.root']

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.04)
c.SetBottomMargin(0.12)


for plot in PLOTS:

    totalMC=None
    for url,norm in mc:

        if '1l4j2b' in plot   : norm *= 0.65**2
        if '1l4j1b1q' in plot : norm *= 2*0.65*(1-0.65)
        if '1l4j2q' in plot   : norm *= (1-0.65)**2
        
        if '_b_' in plot : norm *= 2

        fIn=ROOT.TFile.Open(url)
        h=fIn.Get(plot)
        h.Scale(norm/h.Integral(0,h.GetNbinsX()+1))
        if totalMC : 
            totalMC.Add(h)
        else :
            totalMC=h.Clone(plot+'_mc')
            totalMC.SetDirectory(0)
            totalMC.SetFillStyle(1001)
            ci=ROOT.TColor.GetColor('#91bfdb')
            totalMC.SetFillColor(ci)
            totalMC.SetLineColor(ci)
        fIn.Close()

    totalData=None
    for url in data:
        fIn=ROOT.TFile.Open(url)
        h=fIn.Get(plot)
        if totalData :
            totalData.Add(h)
        else :
            totalData=h.Clone(plot+'_data')
            totalData.SetDirectory(0)
    fIn.Close()

    c.Clear()
    totalMC.Draw('hist')
    totalData.Draw('epsame')
    totalMC.GetYaxis().SetRangeUser(0,totalData.GetMaximum()*1.3)
    leg = ROOT.TLegend(0.14,0.88,0.5,0.75)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.AddEntry(totalData,'data','ep')
    leg.AddEntry(totalMC,'PYQUEN t#bar{t}','f')
    leg.Draw()
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.04)
    label.DrawLatex(0.15,0.9,'#bf{CMS} #it{preliminary}')
    label.DrawLatex(0.58,0.96,'#scale[0.8]{180 nb^{-1} (pPb at #sqrt{s}=8.16 TeV)}')
    c.cd()
    c.Modified()
    c.Update()
    c.SaveAs('%s.pdf'%plot)

