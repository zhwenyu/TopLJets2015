#!/usr/bin/env python

import ROOT

def getCompHistos(COMP,sumCats=True):

    histos=[]
    for url,hname in COMP:
        fIn=ROOT.TFile.Open(url)
        h=None
        for ch in ['EM','MM','EE'] if sumCats else ['']:
            for b in ['1b','2b'] if sumCats else ['']:
                for cat in ['lowpt_','highpt_'] if sumCats else ['']:
                    ih=fIn.Get(ch+b+cat+hname)
                    if h:
                        h.Add( ih )
                    else:
                        h=ih.Clone()
                        h.SetDirectory(0)
        histos.append(h)
        histos[-1].Scale(1./histos[-1].Integral(0,histos[-1].GetNbinsX()+1))
        if len(histos)>1 : 
            histos[-1].Divide(histos[0])
        fIn.Close()
    
    return histos[-1]


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
#ROOT.gROOT.SetBatch(True)


name='mlbweighting_closure'
sumCats=True
COMPLIST={
    '#frac{1}{5}#Gamma_{SM}' : [('analysis/MC13TeV_TTJets_widthx0p2.root','incmlb_w100'), ('analysis/MC13TeV_TTJets.root','incmlb_w20')],
    '#frac{1}{2}#Gamma_{SM}' : [('analysis/MC13TeV_TTJets_widthx0p5.root','incmlb_w100'), ('analysis/MC13TeV_TTJets.root','incmlb_w50')],
    '4#Gamma_{SM}'           : [('analysis/MC13TeV_TTJets_widthx4.root','incmlb_w100'),   ('analysis/MC13TeV_TTJets.root','incmlb_w400')],
    }

name='tmassweighting_closure'
sumCats=False
COMPLIST={
    '#frac{1}{5}#Gamma_{SM}' : [('analysis/MC13TeV_TTJets_widthx0p2.root','tmass_w100'), ('analysis/MC13TeV_TTJets.root','tmass_w20')],
    '#frac{1}{2}#Gamma_{SM}' : [('analysis/MC13TeV_TTJets_widthx0p5.root','tmass_w100'), ('analysis/MC13TeV_TTJets.root','tmass_w50')],
    '4#Gamma_{SM}'           : [('analysis/MC13TeV_TTJets_widthx4.root',  'tmass_w100'), ('analysis/MC13TeV_TTJets.root','tmass_w400')],
    }


histos=[]
for title in COMPLIST:
    histos.append( getCompHistos(COMPLIST[title],sumCats=sumCats) )
    histos[-1].SetTitle(title)
print histos

c= ROOT.TCanvas("c","c",500,500)
c.SetTopMargin(0.06)
c.SetRightMargin(0.03)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.1)

leg=ROOT.TLegend(0.15,0.93,0.6,0.69)
leg.SetTextFont(42)
leg.SetTextSize(0.03)
leg.SetFillStyle(0)
leg.SetFillColor(0)
leg.SetBorderSize(0)

for i in xrange(0,len(histos)):
    histos[i].SetFillStyle(3003+i)
    histos[i].SetFillColor(i*10+20)
    histos[i].SetMarkerColor(i*10+20)
    histos[i].SetMarkerStyle(24+i)
    histos[i].Fit('pol1','+','same')
    pol1=histos[i].GetFunction('pol1')
    pol1.SetLineColor(i*10+20)
    histos[i].Draw('pe3' if i==0 else 'pe3same')
    leg.AddEntry(histos[i],histos[i].GetTitle(),'f')
    histos[i].GetYaxis().SetRangeUser(0.9,1.1)
    histos[i].GetYaxis().SetTitleOffset(1.4)
    histos[i].GetYaxis().SetTitle('Reweighted SM to simulated ratio')
    if 'tmass' in name:
        histos[i].GetXaxis().SetRangeUser(168,177)
leg.Draw()

tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
tex.DrawLatex(0.85,0.96,'13 TeV')

c.Modified()
c.Update()
for ext in ['png','pdf']:
    c.SaveAs('%s.%s'%(name,ext))
