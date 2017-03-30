#!/usr/bin/env/python

import sys
import os
import optparse
import ROOT
from UEAnalysisHandler import *

"""
"""
def getPurStab(h):
    
    for xbin in range(1, h.GetNbinsX()+1):
        for ybin in range(1, h.GetNbinsY()+1,2):
            gensum += h.GetBinContent(g, r)
        gensums.append(gensum)
    
    for r in range(1, h.GetNbinsY()+1,2):
        recosum = 0
        for g in range(0, h.GetNbinsX()+1):
            recosum += h.GetBinContent(g, r)
            recosum += h.GetBinContent(g, r+1)
        recosums.append(recosum)

    
    purGr=ROOT.TGraph()
    stabGr=ROOT.TGraph()
    purGr.SetMarkerStyle(20)
    purGr.SetTitle('purity')
    stabGr.SetMarkerStyle(24)
    stabGr.SetTitle('stability')
    for i in range(len(gensums)):
        purity    = -1
        stability = -1
        if recosums[i] > 0:
            purity = diagonal[i]/recosums[i]
        purGr.SetPoint(purGr.GetN(),i+1,purity)
        if gensums[i] > 0:
            stability = diagonal[i]/gensums[i]
        stabGr.SetPoint(stabGr.GetN(),i+1,stability)
    return purGr,stabGr
    
    

"""
"""
def showMatrices(opt):

    #prepare canvas
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
    ROOT.gStyle.SetPaintTextFormat("4.0f");
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.15)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.1)

    #open ROOT file with the definitions
    ueHandler=UEAnalysisHandler(opt.cfg)
    fIn=ROOT.TFile.Open(opt.input)

    #build the fully combined migration matrices
    for key in ueHandler.histos:
        if len(key)!=2 : continue
        var,a=key
        try:
            stit=VARS[a][0] if a!='inc' else 'inclusive'
        except:
            continue
        otit=VARS[var][0]
        if stit!='inclusive' : continue

        #get histogram and normalize it in gen slices
        hname=ueHandler.histos[key].GetName()
        h2d=fIn.Get(hname)
        for xbin in xrange(1,h2d.GetNbinsX()+1):
            tmp=h2d.ProjectionY('tmp',xbin,xbin)
            total=tmp.Integral(1,tmp.GetNbinsX())
            tmp.Delete()
            if total==0 : continue
            for ybin in xrange(1,h2d.GetNbinsY()+1):
                val=h2d.GetBinContent(xbin,ybin)
                unc=h2d.GetBinError(xbin,ybin)
                h2d.SetBinContent(xbin,ybin,100.*val/total)
                h2d.SetBinError(xbin,ybin,100.*unc/total)

                
        #display the matrix
        c.Clear()
        h2d.Draw('colz text')
        h2d.GetZaxis().SetLabelSize(0.03)
        h2d.GetYaxis().SetLabelSize(0.03)
        h2d.GetXaxis().SetLabelSize(0.03)
        h2d.GetXaxis().SetTitle('Generator level bin')
        h2d.GetYaxis().SetTitle('Reconstruction level bin')
        h2d.GetZaxis().SetRangeUser(0.,100.)

        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.04)
        tex.SetNDC()
        tex.DrawLatex(0.12,0.90,'#bf{CMS} #it{simulation preliminary}')
        tex.DrawLatex(0.12,0.85,otit)
        tex.DrawLatex(0.68,0.96,'#scale[0.8]{36.7 fb^{-1} (13 TeV)}')

        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('%s/%s_migration.%s'%(opt.out,var,ext))


"""
"""
def main():

    #config
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',    dest='input',  help='input',   default='./UEanalysis/analysis_0_0/MC13TeV_TTJets.root', type='string')
    parser.add_option('-c', '--cfg',   dest='cfg',    help='cfg',     default='./UEanalysis/analysiscfg.pck',    type='string')
    parser.add_option('-o', '--out',   dest='out',    help='output',  default='./UEanalysis',                    type='string') 
    (opt, args) = parser.parse_args()
    showMatrices(opt)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
