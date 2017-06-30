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
    gensumswithuf = []

    for g in range(1, h.GetNbinsX()+1):
        gensum = 0
        indices.append(g)
        diagonal.append(h.GetBinContent(g, g))
        for r in range(1, h.GetNbinsY()+1):
            gensum += h.GetBinContent(g, r)
        gensums.append(gensum)
        gensumswithuf.append(gensum+h.GetBinContent(g,0))

    for r in range(1, h.GetNbinsY()+1):
        recosum = 0
        for g in range(0, h.GetNbinsX()+1):
            recosum += h.GetBinContent(g, r)
        recosums.append(recosum)

    
    purGr=ROOT.TGraph()
    stabGr=ROOT.TGraph()
    effGr=ROOT.TGraph()

    for i in range(len(gensums)):
        purity    = -1
        stability = -1
        eff       = -1
        if recosums[i] > 0:
            purity = diagonal[i]/recosums[i]
        purGr.SetPoint(purGr.GetN(),i+1,purity)
        if gensums[i] > 0:
            stability = diagonal[i]/gensums[i]
        stabGr.SetPoint(stabGr.GetN(),i+1,stability)
        if gensumswithuf[i]>0:
            eff=gensums[i]/gensumswithuf[i]
        effGr.SetPoint(effGr.GetN(),i+1,eff)
    return purGr,stabGr,effGr
    
    

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

    c2=ROOT.TCanvas('c2','c2',500,500)
    c2.SetTopMargin(0.05)
    c2.SetRightMargin(0.05)
    c2.SetBottomMargin(0.1)
    c2.SetLeftMargin(0.1)


    #open ROOT file with the definitions
    ueHandler=UEAnalysisHandler(opt.cfg)
    fIn=ROOT.TFile.Open(opt.input)

    #build the fully combined migration matrices
    conditionNumbers=[]
    for key in ueHandler.histos:

        obs,sliceVar,sliceIdx,sel,systVar,distType=key
        title=VARS[obs][0]
        sliceVarTitle=VARS[sliceVar][0] if not sliceVar is None else ''
        outName='_'.join(map(str,key))

        if distType!='mig' : continue
        if systVar!=0 : continue
        #if sel != 'inc' : continue

        #get histogram  
        hname=ueHandler.histos[key].GetName()
        h2d=fIn.Get(hname)

        #this has to be done before normalizing per column otherwise efficiency 
        #and SVD decomposition are meaningless
        h2dForTest=h2d.Clone('h2dfortest')
        h2dForTest.RebinY()
        purGr,stabGr,effGr=getPurStab(h2dForTest)

        try:
            nx,ny=h2dForTest.GetNbinsX(),h2dForTest.GetNbinsY()
            matrix=ROOT.TMatrixD(ny+1,nx)
            for xbin in xrange(1,nx+1):
                for ybin in xrange(0,ny+1):
                    matrix[ybin][xbin-1]=h2dForTest.GetBinContent(xbin,ybin)
            svd=ROOT.TDecompSVD(matrix)
            sig=svd.GetSig()
            maxSig,minSig=sig.Max(),sig.Min()
            condK=-1 if minSig==0 else maxSig/max(0,minSig)
            conditionNumbers.append( (hname,maxSig,minSig,condK) )
        except:
            print 'unable to run SVDdecomposition for',hname

        #normalize it in gen slices
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
        c.cd()
        c.Clear()
        if sliceVar:
            h2d.Draw('colz')
        else:
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
        if sliceVar:
            tex.DrawLatex(0.12,0.85,title+', as function of '+sliceVarTitle)
        else:
            tex.DrawLatex(0.12,0.85,title)
        tex.DrawLatex(0.68,0.96,'#scale[0.8]{35.9 fb^{-1} (13 TeV)}')

        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('%s/%s_migration.%s'%(opt.out,outName,ext))

        # purity efficiency and all that
        c2.cd()
        c2.Clear()
        stabGr.SetLineWidth(2)
        stabGr.SetTitle('stability')
        stabGr.SetName('stab')
        stabGr.Draw('al')
        stabGr.GetYaxis().SetRangeUser(0,1)
        stabGr.GetXaxis().SetTitle('Generator level bin')
        stabGr.GetYaxis().SetTitle('Purity,stability or efficiency')
        purGr.SetLineWidth(2)
        purGr.SetLineColor(ROOT.kRed)
        purGr.SetTitle('purity')
        purGr.SetName('pur')
        purGr.Draw('l')
        effGr.SetLineWidth(2)
        effGr.SetLineColor(ROOT.kGray)
        effGr.SetLineStyle(2)
        effGr.SetTitle('efficiency')
        effGr.SetName('eff')
        effGr.Draw('l')

        leg=ROOT.TLegend(0.12,0.4,0.4,0.2)
        leg.SetTextFont(42)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.AddEntry(stabGr, 'stability',  'l')
        leg.AddEntry(purGr,  'purity',     'l')
        leg.AddEntry(effGr,  'efficiency', 'l')
        leg.Draw()

        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.04)
        tex.SetNDC()
        tex.DrawLatex(0.12,0.90,'#bf{CMS} #it{simulation preliminary}')
        if sliceVar:
            tex.DrawLatex(0.12,0.85,title+', as function of '+sliceVarTitle)
        else:
            tex.DrawLatex(0.12,0.85,title)
        tex.DrawLatex(0.68,0.96,'#scale[0.8]{35.9 fb^{-1} (13 TeV)}')

        c2.cd()
        c2.Modified()
        c2.Update()
        for ext in ['png','pdf']:
            c2.SaveAs('%s/%s_purstabeff.%s'%(opt.out,outName,ext))

        #raw_input()
    return conditionNumbers



"""
"""
def main():

    #config
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',    dest='input',  help='input',   default='./UEanalysis/analysis/MC13TeV_TTJets.root', type='string')
    parser.add_option('-c', '--cfg',   dest='cfg',    help='cfg',     default='./UEanalysis/analysiscfg.pck',    type='string')
    parser.add_option('-o', '--out',   dest='out',    help='output',  default='./UEanalysis',                    type='string') 
    (opt, args) = parser.parse_args()
    conditionNumbers=showMatrices(opt)

    with open('%s/conditionnumbers.dat'%opt.out,'w') as f:
        for k in conditionNumbers: 
            if k[3]<0: continue
            f.write('%25s & %3.5f & %3.5f & %3.2f\n'%k)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
