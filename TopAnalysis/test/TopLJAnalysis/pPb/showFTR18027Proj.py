#!/usr/bin/env python

import ROOT
C90=1.64485
onlyStatFromSplot=True


def getMCFMPreds(migF):

    mcfmPreds={}
    for d in ['pt_l','y_l']:
        mig=getSummed(migF,d+'_mig')
        f=ROOT.TFile.Open('{0}_CT14_nCTEQ15_PDFunc.root'.format(d))
        p=ROOT.TGraph(f.Get(d).At(0))

        genX=mig.ProjectionY('{0}_gen'.format(d))
        genX.Reset('ICE')
        smearedX=mig.ProjectionX('{0}_smaeared'.format(d))
        smearedX.Reset('ICE')

        for xbin in xrange(1,genX.GetNbinsX()+1):
            xcen=genX.GetXaxis().GetBinCenter(xbin)
            val=p.Eval(xcen)
            #if d=='y_l':
            #    val=p.Eval( -1*xcen ) 
            genX.Fill(xcen,val)
            
        for xbin in xrange(1,genX.GetNbinsX()+1):
            
            wgtSum=[]
            for ybin in xrange(1,smearedX.GetNbinsX()+1):
                wgtSum.append(mig.GetBinContent(xbin,ybin))
            totalWgts=sum(wgtSum)
            
            for ybin in xrange(1,smearedX.GetNbinsX()+1):                         
                smearedX.Fill(smearedX.GetXaxis().GetBinCenter(ybin),genX.GetBinContent(xbin)*wgtSum[ybin-1]/totalWgts)
                
        for xbin in xrange(1,genX.GetNbinsX()+1):
            genX.SetBinContent(xbin,genX.GetBinContent(xbin)/genX.GetXaxis().GetBinWidth(xbin))
        for xbin in xrange(1,smearedX.GetNbinsX()+1):
            smearedX.SetBinContent(xbin,smearedX.GetBinContent(xbin)/smearedX.GetXaxis().GetBinWidth(xbin))

        mcfmPreds[d]=smearedX.Clone(d+'_mcfm')
        mcfmPreds[d].SetDirectory(0)
        mcfmPreds[d].SetTitle('MCFM CT14+nCTEQ15')
        mcfmPreds[d].SetLineColor(ROOT.kRed)
        mcfmPreds[d].SetLineWidth(2)
    return mcfmPreds


def cmsHeader(doLumi=False):
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{Simulation}')
    if doLumi:
        tex.SetTextAlign(32)
        tex.DrawLatex(0.95,0.97,'#scale[0.8]{pPb (2pb^{-1}, #sqrt{s}=8.16 TeV)}')


def getPurStab(h):
    """
    get purity, stability and efficiency of a migration matrix
    """
    indices  = []
    diagonal = []
    gensums  = []
    recosums = []
    gensumswithuf = []

    xcen=[]
    for g in range(1, h.GetNbinsX()+1):
        gensum = 0
        indices.append(g)
        diagonal.append(h.GetBinContent(g, g))
        for r in range(1, h.GetNbinsY()+1):
            gensum += h.GetBinContent(g, r)
        gensums.append(gensum)
        gensumswithuf.append(gensum+h.GetBinContent(g,0))
        xcen.append( h.GetXaxis().GetBinCenter(g) )

    for r in range(1, h.GetNbinsY()+1):
        recosum = 0
        for g in range(0, h.GetNbinsX()+1):
            recosum += h.GetBinContent(g, r)
        recosums.append(recosum)

    
    purGr=ROOT.TGraph()
    purGr.SetLineWidth(2)
    purGr.SetName('pur')
    stabGr=purGr.Clone('stab')
    stabGr.SetLineColor(ROOT.kGreen+3)

    for i in range(len(gensums)):
        purity    = -1
        stability = -1
        if recosums[i] > 0:
            purity = diagonal[i]/recosums[i]
        purGr.SetPoint(purGr.GetN(),xcen[i],purity)
        if gensums[i] > 0:
            stability = diagonal[i]/gensums[i]
        stabGr.SetPoint(stabGr.GetN(),xcen[i],stability)
        
    return purGr,stabGr

def normalized(h):
    for xbin in xrange(1,h.GetNbinsX()+1):
        wid=h.GetXaxis().GetBinWidth(xbin)
        val,valUnc=h.GetBinContent(xbin),h.GetBinError(xbin)
        h.SetBinContent(xbin,val/wid)
        h.SetBinError(xbin,valUnc/wid)


def getSummed(inF,hname,divideXwid=False):
    histTotal=None
    for cat in ['1l4j1b1q','1l4j2b']:
        for ch in ['e','mu']:
            h=inF.Get('{0}_{1}{2}'.format(hname,ch,cat))
            if not histTotal:
                histTotal=h.Clone(hname+'_total')
                histTotal.Reset('ICE')
                histTotal.SetDirectory(0)
            histTotal.Add(h)            
    if divideXwid : normalized(histTotal)
    return histTotal

inF=ROOT.TFile.Open('plots/MC8.16TeV_TTbar_pPb_Pohweg/controlplots.root')
splots=ROOT.TFile.Open('splots.root')

mcfmPreds=getMCFMPreds(inF)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
#ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
ROOT.gStyle.SetPaintTextFormat("4.0f")
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetRightMargin(0.03)
c.SetBottomMargin(0.1)
c.SetLeftMargin(0.12)
for hname in ['pt_l','y_l','pt_l_mig','y_l_mig']:
    if 'mig' in hname:

        #show purity, stability etc
        h=getSummed(inF,hname)
        purGr,stabGr=getPurStab(h)
        stabGr.SetLineWidth(2)
        stabGr.SetTitle('stability')
        stabGr.SetName('stab')
        stabGr.Draw('al')
        stabGr.GetYaxis().SetRangeUser(0.8,1)
        stabGr.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
        stabGr.GetYaxis().SetTitle('Purity or stability')
        stabGr.GetYaxis().SetTitleOffset(1.4)
        purGr.SetLineWidth(2)
        purGr.SetLineColor(ROOT.kRed)
        purGr.SetTitle('purity')
        purGr.SetName('pur')
        purGr.Draw('l')

        leg=ROOT.TLegend(0.7,0.94,0.9,0.84)
        leg.SetTextFont(42)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.AddEntry(stabGr, 'stability',  'l')
        leg.AddEntry(purGr,  'purity',     'l')
        leg.Draw()

        cmsHeader()

        c.cd()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('{0}_purstab.{1}'.format(hname,ext))
  
        #normalize migration matrix per gen slices
        for xbin in xrange(1,h.GetNbinsX()+1):
            tmp=h.ProjectionY('tmp',xbin,xbin)
            total=tmp.Integral(1,tmp.GetNbinsX())
            tmp.Delete()
            if total==0 : continue
            for ybin in xrange(1,h.GetNbinsY()+1):
                val=h.GetBinContent(xbin,ybin)
                unc=h.GetBinError(xbin,ybin)
                wid=1.0 #h.GetYaxis().GetBinWidth(ybin)*h.GetXaxis().GetBinWidth(xbin)
                h.SetBinContent(xbin,ybin,100.*val/(total*wid))
                h.SetBinError(xbin,ybin,100.*unc/total)
        h.Draw('col text')
        h.GetYaxis().SetTitleOffset(1.3)
        h.GetXaxis().SetTitle('Generated %s'%h.GetXaxis().GetTitle())
        h.GetYaxis().SetTitle('Reconstructed %s'%h.GetYaxis().GetTitle())
        cmsHeader()
        c.cd()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('{0}.{1}'.format(hname,ext))

    else:
        ytitle='p_{T} [GeV^{-1}]' if hname=='pt_l' else 'y'

        h=getSummed(inF,hname,True)
        sh=splots.Get('splot_{0}__{0}'.format('lpt' if 'pt' in hname else 'ly'))
        sh.SetLineWidth(3)
        sf=sh.Integral()/h.Integral()
        h.Scale(sf)
        #sf=59.*174.5/200000.
        #sf=59.*2000./200000.
        #sf=1./h.Integral(0,h.GetNbinsX()+1)
        #h.Scale(sf)
        hpdf=getSummed(inF,hname+'_pdf')
        frame=None
        projh=[]
        ratioh=[]
        for ybin in xrange(1,hpdf.GetNbinsY()+1):
            projh.append(hpdf.ProjectionX('{0}_p{1}'.format(hname,ybin),ybin,ybin))
            projh[-1].SetDirectory(0)
            projh[-1].Scale(sf)
            normalized(projh[-1])
            if not frame:
                frame=projh[-1].Clone('frame')                
                frame.Reset('ICE')                
                frame.Draw()
                frame.GetYaxis().SetTitleOffset(1.5)
                frame.GetYaxis().SetTitle('dN/d%s'%ytitle)
                xtit=h.GetXaxis().GetTitle()
                xtit=xtit[0].lower()+xtit[1:]
                frame.GetXaxis().SetTitle('Lepton ' + xtit)

        #PDF variations
        varMax=hpdf.Clone('varMax')
        varMin=hpdf.Clone('varMin')
        varMax.Reset('ICE')
        varMin.Reset('ICE')
        for i in xrange(0,len(projh),2):
            up=projh[i]
            dn=projh[i+1]
            for xbin in xrange(1,up.GetNbinsX()+1):
                valup=projh[i].GetBinContent(xbin)
                valdn=projh[i+1].GetBinContent(xbin)
                valcen=h.GetBinContent(xbin)
                varMax.SetBinContent(xbin,varMax.GetBinContent(xbin)+max(valup-valcen,valdn-valcen)**2)
                varMin.SetBinContent(xbin,varMin.GetBinContent(xbin)+max(valup-valcen,valdn-valcen)**2)
        pdfEnv=ROOT.TGraphAsymmErrors()
        pdfRelEnv=ROOT.TGraphAsymmErrors()
        for xbin in xrange(1,h.GetNbinsX()+1):
            cenVal=h.GetBinContent(xbin)

            #normalize to predicted
            if onlyStatFromSplot:
                sh.SetBinError(xbin,sh.GetBinError(xbin)*cenVal/sh.GetBinContent(xbin))
                sh.SetBinContent(xbin,cenVal)

            pdfEnv.SetPoint(xbin-1,h.GetBinCenter(xbin),cenVal)
            pdfEnv.SetPointError(xbin-1,
                                 0.5*h.GetBinWidth(xbin),
                                 0.5*h.GetBinWidth(xbin),
                                 ROOT.TMath.Sqrt(varMin.GetBinContent(xbin))/C90,
                                 ROOT.TMath.Sqrt(varMax.GetBinContent(xbin))/C90)
            pdfRelEnv.SetPoint(xbin-1,h.GetBinCenter(xbin),1)
            pdfRelEnv.SetPointError(xbin-1,
                                    0.5*h.GetBinWidth(xbin),
                                    0.5*h.GetBinWidth(xbin),
                                    ROOT.TMath.Sqrt(varMin.GetBinContent(xbin))/(C90*cenVal),
                                    ROOT.TMath.Sqrt(varMax.GetBinContent(xbin))/(C90*cenVal))
        for g in [pdfEnv,pdfRelEnv]:
            g.SetFillStyle(3001)
            g.SetFillColor(ROOT.kAzure+7)
            g.SetMarkerStyle(1)
            g.SetMarkerColor(ROOT.kAzure+7)
            g.SetLineColor(1)
        pdfEnv.Draw('2')

        for i in xrange(0,len(projh)):
            projh[i].SetLineWidth(2)
            projh[i].SetLineColor(ROOT.kTeal-8)
            #projh[i].Draw('histsame')
            projh[i].GetYaxis().SetTitle(h.GetYaxis().GetTitle())
            projh[i].GetXaxis().SetTitle(h.GetXaxis().GetTitle())
            ratioh.append(projh[i].Clone('{0}_p{1}_ratio'.format(hname,ybin)))
            ratioh[i].Divide(h)                  

        if hname in mcfmPreds:
            mcfmPreds[hname].Scale(h.Integral()/mcfmPreds[hname].Integral())
            mcfmPreds[hname].Draw('histsame')
        
        h.Draw('histsame')
        ROOT.gStyle.SetErrorX(0.)
        sh.Draw('e1same')
        sh.SetMarkerStyle(20)

        leg=ROOT.TLegend(0.5,0.94,0.9,0.8)
        leg.SetTextFont(42)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.AddEntry(sh,       'Pseudo-data', 'ep')
        leg.AddEntry(pdfEnv,   'pPb#rightarrow t#bar{t} Powheg (EPPS16)',  'lf')
        #leg.AddEntry(projh[0], 'EPPS16 error sets',     'l')
        if hname in mcfmPreds:
            leg.AddEntry(mcfmPreds[hname],mcfmPreds[hname].GetTitle(),'l')

        leg.Draw()
        frame.GetYaxis().SetRangeUser(0,h.GetMaximum()*1.5)
        cmsHeader(True)
        c.RedrawAxis()
        c.cd() 
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('{0}{1}.{2}'.format(hname,'realstat' if onlyStatFromSplot else '',ext))

        fOut=ROOT.TFile.Open(hname+'.root','RECREATE')
        h.Write('central')
        sh.Write('splot')
        for i in xrange(0,len(projh)):
            projh[i].Write('npdf_%d'%i)
        fOut.Close()

        frame.GetYaxis().SetTitle('Relative uncertainty')
        frame.GetYaxis().SetRangeUser(0.8,1.2)
        frame.Draw()
        pdfRelEnv.Draw('2')
        #for r in ratioh:             
        #    r.Draw('histsame')
        for xbin in xrange(1,sh.GetNbinsX()+1):
            val=sh.GetBinContent(xbin)
            valUnc=sh.GetBinError(xbin)
            if val==0: continue
            sh.SetBinContent(xbin,1)
            sh.SetBinError(xbin,valUnc/val)

        sh.SetLineColor(1)
        sh.SetMarkerStyle(20)
        sh.SetMarkerColor(1)
        sh.Draw('e1same')
        #h.Divide(h)
        #h.Draw('e1same')
        leg=ROOT.TLegend(0.12,0.94,0.95,0.9)
        leg.SetTextFont(42)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.SetNColumns(3)
        leg.AddEntry(sh,        'Pseudo-data', 'ep')
        leg.AddEntry(pdfRelEnv, 'pPb#rightarrow t#bar{t} Powheg (EPPS16)',  'lf')
        #leg.AddEntry(ratioh[0], 'EPPS16 error sets',     'l')
        leg.Draw()

        cmsHeader(True)
        c.RedrawAxis()
        c.cd()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('{0}{1}_pdfvars.{2}'.format(hname,'realstat' if onlyStatFromSplot else '',ext))        
