#!/usr/bin/env python

import ROOT
import sys

SYSTLIST={
    'LES':['ees','mes'],
    'JER':['jer'],
    'Selection':['sel_E','sel_M','trig_EM','trig_MM','trig_EE'],
    'b-tag':['ltag','btag'],
    'b-fragmentation':['bfrag','semilep'],
    'Pileup':['pu'],
    'Top p_{T}':['tttoppt'],
    'PS scales':['ISR','FSR'],
    'ME QCD scale':['ttMEqcdscale'],
    'JES':['jes%d'%i for i in xrange(0,29)],
    'PDF':['ttPDF'],
    'm_{top}':['mtop'],
    'np QCD':['UE','CR'],
    'hdamp':['hdamp']
}

COLORS=[ROOT.kBlack, ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kMagenta-9,ROOT.kRed+1,ROOT.kAzure+7, ROOT.kBlue-7]
STYLES=[1,2,3,7,9]


def getVariationRatios(url,cat,systs,proc='tbartw400'):
    """computes the ratios for each shape variation"""
    fIn=ROOT.TFile.Open(url)

    central=fIn.Get(cat+'_incmlb/'+proc)

    histos=[]
    for s in systs:

        color=COLORS[ len(histos)%len(COLORS) ]
        style=STYLES[ len(histos)/len(COLORS) ]

        try:
            upVar=fIn.Get('%s_incmlb_%sUp/%s'%(cat,s,proc))
            upVar.Divide(central)
            upVar.SetName('%sUp'%s)
            upVar.SetLineColor(color)
            upVar.SetMarkerColor(color)
            upVar.SetFillStyle(0)
            upVar.SetLineStyle(style)
            upVar.SetTitle(s)
            upVar.SetDirectory(0)
        except:
            upVar=None
            pass

        try:
            downVar=fIn.Get('%s_incmlb_%sDown/%s'%(cat,s,proc))
            downVar.Divide(central)
            downVar.SetName('%sDown'%s)
            downVar.SetDirectory(0)
            downVar.SetLineColor(color)
            downVar.SetMarkerColor(color)
            downVar.SetFillStyle(0)
            downVar.SetLineStyle(style)
            downVar.SetTitle(s)
        except:
            downVar=None
            pass

        histos.append( (upVar,downVar) )
    
    fIn.Close()
    return histos

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

c= ROOT.TCanvas("c","c",500,500)
c.SetTopMargin(0.06)
c.SetRightMargin(0.03)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.1)


url=sys.argv[1]
cat=sys.argv[2]

catTitle=''
if 'EM' in cat: catTitle='e#mu'
if 'EE' in cat: catTitle='ee'
if 'MM' in cat: catTitle='#mu#mu'
if '2b' in cat: catTitle+=',#geq2b'
if '1b' in cat: catTitle+=',=1b'
if 'highpt' in cat: catTitle+=',p_{T}>100 GeV'
if 'lowpt' in cat: catTitle+=',p_{T}<100 GeV'


for s in SYSTLIST:
    histos=getVariationRatios(url=url,cat=cat,systs=SYSTLIST[s])

    c.Clear()
    leg=ROOT.TLegend(0.15,0.86,0.95,0.64)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)

    filled=False
    for i in xrange(0,len(histos)):
        for h in histos[i]:
            if h is None: continue
            h.Draw('histsame' if filled else 'hist')
            filled=True
            if 'Up' in h.GetName(): leg.AddEntry(h,h.GetTitle(),'l')
            h.GetXaxis().SetTitle('Mass(lepton,b) [GeV]')
            h.GetXaxis().SetRangeUser(20,200)
            h.GetYaxis().SetTitleOffset(1.4)
            h.GetYaxis().SetTitle('Ratio to nominal')
            h.GetYaxis().SetRangeUser(0.85,1.15)
    if len(SYSTLIST[s])>4: leg.SetNColumns(4)
    leg.Draw()

    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.DrawLatex(0.85,0.96,'13 TeV')
    tex.DrawLatex(0.15,0.88,'%s #it{%s}'%(catTitle,s))

    c.Modified()
    c.Update()
    name=s
    for t in [' ','{','}','^','-']:name=name.replace(t,'')

    for ext in ['png','pdf']:
        c.SaveAs('%s_%s.%s'%(cat,name,ext))
