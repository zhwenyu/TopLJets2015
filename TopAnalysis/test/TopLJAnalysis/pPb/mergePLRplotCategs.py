#!/usr/bin/env python

import ROOT
import sys
import pickle

lumi=(174.5,8.725)
#from runDataFit import lumi

indchisquares = pickle.load( open( "chisquare_final.pck", "rb" ) )

plots={'pdfcomp1':None,
       'pdfcomp0':None,
       'data':None,
       'totalpdf':None,
       'totalpdfcont':None}

var=sys.argv[1]
varTitle="m_{jj'} [GeV]"
varRange=(25,300)
if 'mthad' in var:
    varTitle='m_{top} [GeV]'
    varRange=(100,400)
if 'mtlep' in var:
    varTitle='m_{t_{lep}} [GeV]'
    varRange=(100,400)

cat=sys.argv[2]
catTitle = 'e^{#pm}/#mu^{#pm} + #geq4j'
if '2q' in cat   : catTitle += ' (=0b)'
if '1b1q' in cat : catTitle += ' (=1b)'
if '2b' in cat  : catTitle += ' (#geq2b)'

noPulls=False
try:
    if sys.argv[3]=='True': noPulls=True
except:
    pass

maxY=0
totalchisq=0
for ch in ['e','mu']:
    fIn=ROOT.TFile.Open('%s_%s%s_final.root'%(var,ch,cat))

    totalchisq+= 0.5*indchisquares[(var,ch+cat)]

    for p in fIn.Get('c').GetPrimitive('p1').GetListOfPrimitives():
        pname=p.GetName()

        if not pname in plots : continue
        print ch,pname,p.Integral()

        if plots[pname] is None :
            plots[pname]=p.Clone()
        else:

            x,y=ROOT.Double(0),ROOT.Double(0)
            if pname!='data':
                
                for n in xrange(0,plots[pname].GetN()):
                    plots[pname].GetPoint(n,x,y)
                    plots[pname].SetPoint(n,x,y+p.Eval(x))
            else:

                xcur,ycur=ROOT.Double(0),ROOT.Double(0)
                for n in xrange(0,p.GetN()):
                    p.GetPoint(n,x,y)
                    
                    for k in xrange(0,plots[pname].GetN()):
                        
                        plots[pname].GetPoint(k,xcur,ycur)
                        if xcur!=x : continue 
                    
                        exhi=p.GetErrorXhigh(n)
                        exlo=p.GetErrorXlow(n)
                        eyhi=p.GetErrorYhigh(n)
                        eylo=p.GetErrorYlow(n)

                        exhicur=plots[pname].GetErrorXhigh(k)
                        exlocur=plots[pname].GetErrorXlow(k)
                        eyhicur=plots[pname].GetErrorYhigh(k)
                        eylocur=plots[pname].GetErrorYlow(k)

                        maxY=max(y+ycur,maxY)
                        plots[pname].SetPoint(k,x,y+ycur)
                        plots[pname].SetPointError(k,
                                                   ROOT.TMath.Sqrt(exlo**2+exlocur**2),
                                                   ROOT.TMath.Sqrt(exhi**2+exhicur**2),
                                                   ROOT.TMath.Sqrt(eylo**2+eylocur**2),
                                                   ROOT.TMath.Sqrt(eyhi**2+eyhicur**2))

    fIn.Close()


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

c=ROOT.TCanvas('c','c',800,800)
c.SetTopMargin(0)
c.SetLeftMargin(0)
c.SetRightMargin(0)
c.SetBottomMargin(0)
p1 = ROOT.TPad('p1','p1',0.0,0.0,1.0,1.0) if noPulls else ROOT.TPad('p1','p1',0.0,0.25,1.0,1.0)
p1.SetRightMargin(0.05)
p1.SetLeftMargin(0.12 if noPulls else 0.15)
p1.SetTopMargin(0.06 if noPulls else 0.08)
p1.SetBottomMargin(0.15 if noPulls else 0.02)
p1.Draw()
p1.cd()
plots['totalpdf'].Draw('alf')
plots['pdfcomp1'].SetPoint(plots['pdfcomp1'].GetN(),plots['totalpdf'].GetXaxis().GetXmax(),0)
plots['pdfcomp0'].SetPoint(plots['pdfcomp0'].GetN(),plots['totalpdf'].GetXaxis().GetXmax(),0)
plots['pdfcomp1'].Draw('lf')
plots['pdfcomp0'].Draw('lf')
for i in xrange(0,plots['data'].GetN()):
    plots['data'].SetPointEXhigh(i,0.)
    plots['data'].SetPointEXlow(i,0.)
plots['data'].Draw('ep')
plots['totalpdf'].GetYaxis().SetRangeUser(0.1,maxY*1.3)
plots['totalpdf'].GetYaxis().SetTitle("Events")
plots['totalpdf'].GetYaxis().SetTitleOffset(1.05 if noPulls else 0.8)
plots['totalpdf'].GetYaxis().SetTitleSize(0.06 if noPulls else 0.07)
plots['totalpdf'].GetYaxis().SetLabelSize(0.05 if noPulls else 0.05)
plots['totalpdf'].GetXaxis().SetTitleSize(0.06 if noPulls else 0.)
plots['totalpdf'].GetXaxis().SetLabelSize(0.05 if noPulls else 0.)
plots['totalpdf'].GetXaxis().SetTitleOffset(0.9 if noPulls else 0.8)
plots['totalpdf'].GetXaxis().SetTitle(varTitle if noPulls else "")
plots['totalpdf'].GetXaxis().SetRangeUser(varRange[0],varRange[1])

label = ROOT.TLatex()
label.SetNDC()
label.SetTextFont(42)
label.SetTextSize(0.05 if noPulls else 0.06)
label.DrawLatex(0.17 if noPulls else 0.18,0.87 if noPulls else 0.85,'#scale[1.2]{#bf{CMS}}') # #it{preliminary}')                                                                     
label.DrawLatex(0.33 if noPulls else 0.37,0.955 if noPulls else 0.95,'pPb (%3.0f nb^{-1}, #sqrt{s_{NN}} = 8.16 TeV)'%lumi[0])
label.DrawLatex(0.50,0.87 if noPulls else 0.86,'#bf{%s}'%catTitle)

pullGr=plots['data'].makeResidHist(plots['totalpdf'],True,True)
ndof=pullGr.GetN()
#if not noPulls : 
#label.DrawLatex(0.50,0.49,'#chi^{2}/dof = %3.1f/%d'%(totalchisq*ndof,ndof))

leg=ROOT.TLegend(0.49,0.85 if noPulls else 0.83,0.8,0.56)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.05 if noPulls else 0.06)
leg.AddEntry(plots['data'],'Data','ep')
leg.AddEntry(plots['totalpdf'],'t#bar{t} correct','f')
leg.AddEntry(plots['pdfcomp1'],'t#bar{t} wrong','f')
leg.AddEntry(plots['pdfcomp0'],'background','f')
leg.Draw()

p1.cd()
p1.RedrawAxis()


c.cd()
if not noPulls:
    p2 = ROOT.TPad('p2','p2',0.0,0.0,1.0,0.25)
    p2.SetBottomMargin(0.45)
    p2.SetRightMargin(0.05)
    p2.SetLeftMargin(0.15)
    p2.SetTopMargin(0.001)
    p2.SetGridy(True)
    p2.Draw()
    p2.cd()
    p2.Clear()
    pullGr.Draw('ap')
    pullGr.GetYaxis().SetRangeUser(-2.4,2.4)
    pullGr.GetYaxis().SetTitle("#frac{Data-Fit}{Unc.}")
    pullGr.GetYaxis().SetTitleSize(0.20)
    pullGr.GetYaxis().SetLabelSize(0.16)
    pullGr.GetYaxis().SetTitleOffset(0.27)
    pullGr.GetXaxis().SetTitleSize(0.20)
    pullGr.GetXaxis().SetLabelSize(0.16)
    pullGr.GetXaxis().SetTitleOffset(0.85)
    pullGr.GetYaxis().SetNdivisions(4)
    pullGr.GetYaxis().SetRangeUser(-3.0,3.0)
    pullGr.GetXaxis().SetTitle(varTitle)
    pullGr.GetXaxis().SetRangeUser(varRange[0],varRange[1])
    p2.cd()
    p2.RedrawAxis()


pullpf='_nopull' if noPulls else ''
c.cd()
c.Modified()
c.Update()
for ext in ['png','pdf','root']:
    c.SaveAs('%s_comb_%s%s.%s'%(var,cat,pullpf,ext))      
p1.cd()
label.DrawLatex(0.17 if noPulls else 0.18,0.8,'#scale[1.0]{#it{preliminary}}')
c.Modified()
c.Update()
for ext in ['png','pdf','root']:
    c.SaveAs('%s_comb_%s_prel%s.%s'%(var,cat,pullpf,ext))   
