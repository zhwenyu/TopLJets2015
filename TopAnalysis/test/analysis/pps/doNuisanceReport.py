#!/usr/bin/env python
import os,sys
import ROOT
from optparse import OptionParser
from collections import defaultdict
import numpy as np
from rounding import *

_skipstat=True

def doNuisanceReport(url,npergroup=20,poi='r'):

    """compare postfit nuisances"""

    fitResults=defaultdict(list)
    poiResults={}
    inF=ROOT.TFile.Open(url)
    for fit in ['b','s']:
        fres=inF.Get('fit_%s'%fit)
        pars=fres.floatParsFinal()
        iter = pars.createIterator()
        var = iter.Next()
        while var :
            name=var.GetName()

            if _skipstat and name.find('prop_bincat')==0:
                var = iter.Next()
                continue
           
            name=name.replace('_bkgShapeSinglediffCat','Single diff. bkg cat')
            name=name.replace('_bkgShapeCat','Bkg. shape cat')
            name=name.replace('mu_outfidsig',"#mu(out-fid.)")
            name=name.replace('mu_bkgCat',"#mu(bkg) cat")
            name=name.replace('sigPzModel','p_{z} signal')
            name=name.replace('eff_','eff.')
            name=name.replace('sigPPSEff','PPS eff.')
            name=name.replace('sigCalib','Time dependency')
            name=name.replace('sigShape','Signal shape')            
            name=name.replace('prop_bincat','stat. cat')
            name=name.replace('_bin',' b')
            
            val=var.getVal()
            ehi=var.getErrorHi()
            elo=var.getErrorLo()
            if name!=poi:
                corr=fres.correlation(var.GetName(),poi) if fit=='s' else 0.
                fitResults[fit].append((name,val,ehi,elo,corr))
            else:
                poiResults[fit]=(val,ehi,elo)

            var = iter.Next()
        
        #sort by constraining power
        fitResults[fit] = sorted(fitResults[fit], key=lambda x: (x[2]*x[2]+x[3]*x[3]))

    #show nuisances
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0)
    c.SetTopMargin(0)
    c.SetRightMargin(0)
    c.SetBottomMargin(0)

    c.cd()
    p1=ROOT.TPad('p1','p1',0,0,0.7,1)
    p1.SetLeftMargin(0.3)
    p1.SetTopMargin(0.1)
    p1.SetRightMargin(0.03)
    p1.SetBottomMargin(1.0)
    p1.SetGridy(True)    
    p1.Draw()

    c.cd()
    p2=ROOT.TPad('p2','p2',0.7,0,1.0,1)
    p2.SetLeftMargin(0.02)
    p2.SetTopMargin(0.1)
    p2.SetRightMargin(0.03)
    p2.SetBottomMargin(1.0)
    p2.SetGridy(True)    
    p2.Draw()
    
    varNames=[x[0] for x in fitResults['s']]
    ngroups=len(varNames)/npergroup
    if ngroups*npergroup<len(varNames) : 
        ngroups+=1

    for ig in range(ngroups):

        first=npergroup*ig
        last=min(npergroup*(ig+1),len(varNames))
        varList=varNames[first:last]
        print(varList)
        #prepare a frame
        npars=len(varList)

        p1.cd()
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
        
        nuisGrs={}
        corrGrs={}
        for iv in range(len(varList)):

            v=varList[iv]

            for fit in fitResults:

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
                    corrGrs[fit]=nuisGrs[fit].Clone(fit+'corr')
                    corrGrs[fit].SetFillStyle(1001)
                    corrGrs[fit].SetFillColor(ROOT.kGray)

                npts=nuisGrs[fit].GetN()
                val,uncLo,uncHi,rho = [x[1:] for x in fitResults[fit] if x[0]==v][0]
                y0=frame.GetYaxis().GetBinCenter(iv+1)
                dy=frame.GetYaxis().GetBinWidth(iv)

                corrGrs[fit].SetPoint(npts,rho,y0)
                corrGrs[fit].SetPointError(npts,rho,0,0.5*dy,0.5*dy)

                if fit=='s': y0-=0.1*dy
                if fit=='b': y0+=0.1*dy
                nuisGrs[fit].SetPoint(npts,val,y0)
                nuisGrs[fit].SetPointError(npts,abs(uncLo),abs(uncHi),0.,0.)



        leg=ROOT.TLegend(0.6,0.94,0.95,0.97)
        leg.SetTextFont(42)
        leg.SetNColumns(2)
        leg.SetTextSize(0.035)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        for fit in nuisGrs:
            nuisGrs[fit].Draw('p')
            leg.AddEntry(nuisGrs[fit],fit,'p')
        leg.Draw()
        p1.RedrawAxis()

        p2.cd()
        frame2=ROOT.TH2F('frame2',';#rho;',1,-0.5,0.5,npars,0,npars)
        frame2.Draw()
        frame2.GetYaxis().SetTitleSize(0)
        frame2.GetYaxis().SetLabelSize(0)
        frame2.GetXaxis().SetTitleSize(0.12)
        frame2.GetXaxis().SetLabelSize(0.09)
        frame2.GetXaxis().SetTitleOffset(0.4)
        frame2.GetXaxis().SetLabelOffset(-0.04)
        frame2.GetXaxis().SetNdivisions(5)
        frame2.GetXaxis().CenterTitle()
        corrGrs['s'].Draw('2')


        txt2=ROOT.TLatex()
        txt2.SetNDC(True)
        txt2.SetTextFont(42)
        txt2.SetTextSize(0.09)
        txt2.SetTextAlign(22)
        poiFit,uncUp,uncDown=poiResults['s']
        txt2.DrawLatex(0.5,0.95,'#hat{%s} = %s'%(poi,toROOTRounded(poiFit,((uncUp,uncDown),))))

        p2.RedrawAxis()

        c.cd()
        c.Modified()
        c.Update()

        for ext in ['png','pdf']:
            c.SaveAs('nuisances_%d.%s'%(ig,ext))

    return

def main():

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    url=sys.argv[1]    
    doNuisanceReport(url=url)


if __name__ == "__main__":
    sys.exit(main())
