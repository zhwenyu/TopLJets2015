#!/usr/bin/env python

import ROOT
import sys
import os

def getEffPlotsFrom(url,tag,outDir):
    fIn=ROOT.TFile.Open(url)
    tpTree=fIn.Get("tpTree")
    effGrs={}
    for key in tpTree.GetListOfKeys():
        name=key.GetName()
        obj=tpTree.Get(name)

        for subkey in obj.GetListOfKeys():
            subname=subkey.GetName()
            subobj=obj.Get(subname)
            try:
                for subsubkey in subobj.GetListOfKeys():
                    subsubname=subsubkey.GetName()
                    cnv=subobj.Get(subsubname)                      
                    if subsubname=='fit_canvas':
                        cnv.Draw()
                        for ext in ['png','pdf']:
                            cnv.SaveAs('%s/%s_%s_%s_%s.%s'%(outDir,tag,name,subname,subsubname,ext))
                    if ('PLOT' in subsubname and subname=='fit_eff_plots'):
                        gr=cnv.GetPrimitive('hxy_fit_eff')
                        try :
                            if gr.GetN()<2 : continue
                            gr.SetTitle(tag)
                            effGrs[(name,subsubname)]=gr
                        except:
                            pass
                        
            except:
                pass
    fIn.Close()
    return effGrs

def getScaleFactor(gr1,gr2):
    sfgr=ROOT.TGraphAsymmErrors()
    x1, y1, x2, y2 = ROOT.Double(0), ROOT.Double(0), ROOT.Double(0), ROOT.Double(0)
    sfmin,sfmax=9999,-9999
    gr1.GetPoint(0,x1,y1)
    xmin=float(x1)-gr1.GetErrorXlow(0)
    gr1.GetPoint(gr1.GetN()-1,x1,y1)
    xmax=float(x1)+gr1.GetErrorXhigh(gr1.GetN()-1)
    for ip in xrange(0,gr1.GetN()):
        gr1.GetPoint(ip,x1,y1)
        gr2.GetPoint(ip,x2,y2)
        exlo=gr2.GetErrorXlow(ip)
        exhi=gr2.GetErrorXhigh(ip)
        ey1=gr1.GetErrorY(ip)
        ey2=gr2.GetErrorY(ip)
        sf=y2/y1
        sferr=sf*ROOT.TMath.Sqrt((ey2/y2)**2+(ey1/y1)**2)
        sfgr.SetPoint(ip,x1,sf)
        sfgr.SetPointError(ip,exlo,exhi,sferr,sferr)
    
        if sferr/sf < 0.1:
            sfmin=ROOT.TMath.Min(sf,sfmin)
            sfmax=ROOT.TMath.Max(sf,sfmax)

    envGr=ROOT.TGraph()
    envGr.SetPoint(0,xmin,sfmin)
    envGr.SetPoint(1,xmin,sfmax)
    envGr.SetPoint(2,xmax,sfmax)
    envGr.SetPoint(3,xmax,sfmin)
    envGr.SetPoint(4,xmin,sfmin)

    finalenvGr=ROOT.TGraph()
    finalenvGr.SetPoint(0,xmin,0.97)
    finalenvGr.SetPoint(1,xmin,1.03)
    finalenvGr.SetPoint(2,xmax,1.03)
    finalenvGr.SetPoint(3,xmax,0.97)
    finalenvGr.SetPoint(4,xmin,0.97)
    return sfgr,envGr,finalenvGr

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

#get efficiency from files
files=sys.argv[1].split(',')
outDir=sys.argv[2]
os.system('mkdir -p %s'%outDir)

effGrs={}
for urltag in files:
    url,tag=urltag.split(':')
    effGrs[tag]=getEffPlotsFrom(url,tag,outDir)

c=ROOT.TCanvas('c','c',500,500)
for item in effGrs['Data']:

    grName=''.join(item)

    grData=effGrs['Data'][item]
    grMC=effGrs['MC'][item]
    sf,envGr,finalenvGr=getScaleFactor(grMC,grData)
    sf.SetLineColor(ROOT.kRed)
    sf.SetMarkerStyle(21)
    sf.SetMarkerColor(ROOT.kRed)
    sf.SetFillStyle(0)
    sf.GetYaxis().SetTitle('Efficiency or scale factor')
    if '_pt' in grName:
        sf.GetXaxis().SetTitle('Transverse momentum [GeV]')
    else:
        sf.GetXaxis().SetTitle('Pseudo-rapidity')
    envGr.SetFillStyle(3244)
    envGr.SetFillColor(ROOT.kRed-9)
    envGr.SetLineColor(ROOT.kRed-9)
    finalenvGr.SetFillStyle(3004)
    finalenvGr.SetFillColor(ROOT.kGray+3)
    finalenvGr.SetLineColor(ROOT.kGray+3)
    grData.SetMarkerStyle(20)
    grData.SetFillStyle(0)
    grMC.SetMarkerStyle(24)
    grMC.SetMarkerColor(ROOT.kGray)
    grMC.SetLineColor(ROOT.kGray)
    grMC.SetFillStyle(0)
    
    c.SetGridx()
    c.SetGridy()

    sf.Draw('ap')
    sf.GetYaxis().SetRangeUser(0.7,1.2)
    grMC.Draw('p')
    grData.Draw('p')    
    sf.SetTitle('Data/MC')
    c.BuildLegend(0.4,0.2,0.8,0.4,'#bf{CMS} #it{preliminary} #sqrt{s_{NN}}=8.16 TeV')
    finalenvGr.Draw('f')
    envGr.Draw('f')
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s/%s.%s'%(outDir,grName,ext)) 
    


