#!/usr/bin/env python
import os,sys
import ROOT
from optparse import OptionParser
from collections import defaultdict
import numpy as np

def doNuisanceReport(url,npergroup=20, blackList=['r']):

    """compare postfit nuisances"""

    fitResults=defaultdict(dict)
    inF=ROOT.TFile.Open(url)
    for fit in ['b','s']:
        fres=inF.Get('fit_%s'%fit)
        pars=fres.floatParsFinal()
        iter = pars.createIterator()
        var = iter.Next()
        while var :
            name=var.GetName()
            if 'prop_' in name:
                name=name.replace('prop_bin','')
                name=name.replace('_bin','')

            if not name in blackList:
                val=var.getVal()
                ehi=var.getErrorHi()
                elo=var.getErrorLo()
                var = iter.Next()
                fitResults[name][fit]=(val,ehi,elo)
            var = iter.Next()

    #show nuisances
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.3)
    c.SetTopMargin(0.1)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(1.0)
    c.SetGridy(True)    
    varNames=fitResults.keys()
    ngroups=len(varNames)/npergroup+1
    for ig in range(ngroups):

        first=npergroup*ig
        last=min(npergroup*(ig+1),len(varNames))
        varList=varNames[first:last]

        #prepare a frame
        npars=len(varList)
        frame=ROOT.TH2F('frame',';#hat{#theta};Nuisance parameter',1,-3,3,npars,0,npars)
        frame.SetDirectory(0)
        for ybin in range(npars):
            frame.GetYaxis().SetBinLabel(ybin+1,'#color[%d]{%s}'%((ybin%2)*10+1,varList[ybin]))
        frame.Draw()
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleOffset(0.8)
        frame.GetXaxis().CenterTitle()
        frame.GetYaxis().SetTitleOffset(3.2)

        #1-sigma
        gr1s=ROOT.TGraph()
        gr1s.SetName('gr1s')
        gr1s.SetMarkerStyle(1)
        gr1s.SetMarkerColor(19)
        gr1s.SetLineColor(19)
        gr1s.SetFillStyle(1001)
        gr1s.SetFillColor(19)
        gr1s.SetPoint(0,-1,0)
        gr1s.SetPoint(1,-1,npars)
        gr1s.SetPoint(2,1,npars)
        gr1s.SetPoint(3,1,0)
        gr1s.SetPoint(4,-1,0)

        #2-sigma
        gr2s=gr1s.Clone('gr2s')
        gr2s.SetMarkerColor(18)
        gr2s.SetLineColor(18)
        gr2s.SetFillStyle(1001)
        gr2s.SetFillColor(18)
        gr2s.SetPoint(0,-2,0)
        gr2s.SetPoint(1,-2,npars)
        gr2s.SetPoint(2,2,npars)
        gr2s.SetPoint(3,2,0)
        gr2s.SetPoint(4,-2,0)

        gr2s.Draw('f')
        gr1s.Draw('f')

        txt=ROOT.TLatex()
        txt.SetTextFont(42)
        txt.SetTextSize(0.03)
        txt.SetTextColor(ROOT.kGray+3)
        txt.SetNDC(True)
        for delta,title in [(0.72,'-1#sigma'),(0.82,'+2#sigma'),(0.5,'-1#sigma'),(0.4,'-2#sigma')]:
            txt.DrawLatex(delta,0.91,title)

        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(42)
        txt.SetTextSize(0.05)
        txt.SetTextAlign(12)
        txt.DrawLatex(0.05,0.955,'#bf{CMS} #it{Preliminary}')
        #txt.SetTextAlign(31)
        #txt.DrawLatex(0.95,0.955,'#scale[0.7]{35.6 fb^{-1} (13 TeV)}')
        
        nuisGrs={}
        for iv in range(len(varList)):

            v=varList[iv]

            for fit in fitResults[v]:

                #start graph if needed
                if not fit in nuisGrs:
                    nuisGrs[fit]=ROOT.TGraphAsymmErrors()
                    nuisGrs[fit].SetTitle(fit)
                    marker=20 if fit=='s' else 24
                    ci=1 if fit=='s' else ROOT.kGreen+2
                    nuisGrs[fit].SetMarkerStyle(marker)
                    nuisGrs[fit].SetMarkerColor(ci)
                    nuisGrs[fit].SetLineColor(ci)
                    nuisGrs[fit].SetLineWidth(2)
                    nuisGrs[fit].SetFillStyle(0)

                npts=nuisGrs[fit].GetN()
                val,uncLo,uncHi = fitResults[v][fit]
                y0=frame.GetYaxis().GetBinCenter(iv+1)
                dy=frame.GetYaxis().GetBinWidth(iv)
                if fit=='s': y0-=0.1*dy
                if fit=='b': y0+=0.1*dy
                nuisGrs[fit].SetPoint(npts,val,y0)
                nuisGrs[fit].SetPointError(npts,abs(uncLo),abs(uncHi),0.,0.)
                
        leg=ROOT.TLegend(0.8,0.9,0.95,0.78)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        for fit in nuisGrs:
            nuisGrs[fit].Draw('p')
            leg.AddEntry(nuisGrs[fit],fit,'p')
        leg.Draw()

        c.RedrawAxis()
        c.Modified()
        c.Update()

        for ext in ['png']: #,'pdf']:
            c.SaveAs('nuisances_%d.%s'%(ig,ext))

    return

def main():

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    doNuisanceReport(url=sys.argv[1])


if __name__ == "__main__":
    sys.exit(main())
