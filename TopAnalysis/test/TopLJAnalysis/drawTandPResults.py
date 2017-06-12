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
c.SetTopMargin(0.05)
c.SetRightMargin(0.05)

for item in effGrs['Data']:

    grName=''.join(item)

    grData=effGrs['Data'][item]
    grData.SetMarkerStyle(20)
    grData.SetFillStyle(0)

    grMCs=[]
    sf=[]
    for tag in effGrs:
        if tag=='Data': continue
        icount=len(sf)

        grMCs.append( effGrs[tag][item] )
        grMCs[-1].SetMarkerStyle(21+icount)
        grMCs[-1].SetMarkerColor(ROOT.kGray+1-icount)
        grMCs[-1].SetLineColor(ROOT.kGray+1-icount)
        grMCs[-1].SetFillStyle(0)
        grMCs[-1].SetTitle(tag)
        sf.append( getScaleFactor(grMCs[-1],grData) )
        sf[-1][0].SetLineColor(ROOT.kRed+icount)
        sf[-1][0].SetMarkerStyle(21+icount)
        sf[-1][0].SetMarkerColor(ROOT.kRed+icount)
        sf[-1][0].SetFillStyle(0)
        sf[-1][0].GetYaxis().SetTitle('Efficiency or scale factor')
        sf[-1][0].SetTitle('Data/%s'%tag)
        if '_pt' in grName:
            sf[-1][0].GetXaxis().SetTitle('Transverse momentum [GeV]')
        else:
            sf[-1][0].GetXaxis().SetTitle('Pseudo-rapidity')
        sf[-1][1].SetFillStyle(3244+1-icount)
        sf[-1][1].SetFillColor(ROOT.kRed-9+icount)
        sf[-1][1].SetLineColor(ROOT.kRed-9+icount)
        sf[-1][2].SetFillStyle(3004+1-icount)
        sf[-1][2].SetFillColor(ROOT.kGray+2-icount)
        sf[-1][2].SetLineColor(ROOT.kGray+2-icount)
        sf[-1][2].SetTitle('Data/%s envelope'%tag)
    
    
    c.SetGridx()
    c.SetGridy()
    drawOpt='ap'
    for i in xrange(0,len(sf)):
        sf[i][0].Draw(drawOpt)        
        sf[i][0].GetYaxis().SetRangeUser(0.6,1.1)
        sf[i][0].GetYaxis().SetTitleOffset(1.1)
        drawOpt='p'
        grMCs[i].Draw(drawOpt)
        #if i==len(sf)-1: 
        #    sf[i][2].Draw('f')
    grData.Draw('p')    

    leg=c.BuildLegend(0.4,0.2,0.8,0.5,'#bf{CMS} #it{preliminary} #sqrt{s_{NN}}=8.16 TeV')
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s/%s.%s'%(outDir,grName,ext)) 
    


