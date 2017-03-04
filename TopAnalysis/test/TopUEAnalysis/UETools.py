#!/usr/bin/env/python

import sys
import os
import optparse
import ROOT
from UEAnalysisHandler import *

"""
"""
def getPurStab(h):
    indices  = []
    diagonal = []
    gensums  = []
    recosums = []
    
    for g in range(1, h.GetNbinsX()+1):
        gensum = 0
        indices.append(g)
        diagonal.append(h.GetBinContent(g, g))
        for r in range(1, h.GetNbinsY()+1):
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
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
    c=ROOT.TCanvas('c','c',1000,1000)
    c.SetTopMargin(0.01)
    c.SetRightMargin(0.1)
    c.SetBottomMargin(0.15)

    
    #open ROOT file with the definitions
    ueHandler=UEAnalysisHandler(opt.cfg)
    fIn=ROOT.TFile.Open(opt.input)

    #build the fully combined migration matrices
    for key in ueHandler.histos:
        if len(key)!=2 : continue
        var,a=key
        stit=VARS[a][0] if a!='inc' else 'inclusive'
        otit=VARS[var][0]
        if stit!='inclusive' : continue

        #get histogram and normalize it in gen slices
        hname=ueHandler.histos[key].GetName()
        h2d=fIn.Get(hname)
        for xbin in xrange(1,h2d.GetNbinsX()+1):
            tmp=h2d.ProjectionY('tmp',xbin,xbin)
            total=tmp.Integral(0,tmp.GetNbinsX()+1)
            tmp.Delete()
            if total==0 : continue
            for ybin in xrange(1,h2d.GetNbinsY()+1):
                val=h2d.GetBinContent(xbin,ybin)
                unc=h2d.GetBinError(xbin,ybin)
                h2d.SetBinContent(xbin,ybin,val/total)
                h2d.SetBinContent(xbin,ybin,unc/total)

        #display the matrix
        c.Clear()
        #c.SetLogz()
        h2d.Draw('colz')
        h2d.GetZaxis().SetLabelSize(0.03)
        h2d.GetYaxis().SetLabelSize(0.03)
        h2d.GetXaxis().SetLabelSize(0.03)
        h2d.GetXaxis().SetTickLength(0)
        h2d.GetYaxis().SetTickLength(0)
        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.025)
        tex.SetNDC()
        tex.DrawLatex(0.01,0.96,'#bf{#splitline{Reco.}{level}}')
        tex.DrawLatex(0.95,0.1,'#bf{#splitline{Gen.}{level}}')
        sep=' ' if stit=='inclusive' else ' vs '
        tex.DrawLatex(0.10,0.03,'#scale[1.2]{#bf{CMS}} #it{simulation preliminary} %s%s%s'%(stit,sep,otit))
        c.Modified()
        c.Update()
        raw_input()
        gr1,gr2=getPurStab(h2d)
        gr1.Draw('ap')
        gr2.Draw('p')
        raw_input()
        #for ext in ['png','pdf']:
        #    c.SaveAs('%s/%s_%s_migration.%s'%(opt.out,s,o,ext))

    #save to ROOT file
#    fOut=ROOT.TFile.Open('%s/UEanalysis.root'%opt.out,'UPDATE')
#    outDir=fOut.Get('final_'+mmDir)
#    try:
#        outDir.cd()
#    except:
#        outDir=fOut.mkdir('final_'+mmDir)
#        outDir.cd()
#    for k in bigMatrix:
#        bigMatrix[k].SetDirectory(outDir)
#        bigMatrix[k].Write(bigMatrix[k].GetName(),ROOT.TObject.kOverwrite)
#    fOut.Close()

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
