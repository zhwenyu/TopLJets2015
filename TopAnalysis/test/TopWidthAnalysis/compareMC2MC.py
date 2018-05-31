#!/usr/bin/env python

import ROOT

def getRatio(num,den,name):
    r=ROOT.TGraphErrors()
    r.SetName(name)
    r.SetTitle(num.GetTitle())
    r.SetMarkerStyle(num.GetMarkerStyle())
    r.SetMarkerColor(num.GetMarkerColor())
    r.SetLineColor(num.GetLineColor())
    r.SetLineWidth(num.GetLineWidth())
    r.SetLineStyle(num.GetLineStyle())

    x,y1,y2=ROOT.Double(0),ROOT.Double(0),ROOT.Double(0)
    for i in xrange(0,num.GetN()):
        num.GetPoint(i,x,y1)
        ey1=num.GetErrorY(i)
        den.GetPoint(i,x,y2)
        ey2=den.GetErrorY(i)
        
        if y1==0: continue
        newy=float(y2)/float(y1)
        if newy==0: continue
        newyUnc=ROOT.TMath.Sqrt((float(y2)*ey1)**2+(float(y1)*ey2)**2)/(float(y1)**2)
                
        npt=r.GetN()
        r.SetPoint(npt,float(x),newy)
        r.SetPointError(npt,0,newyUnc)
    r.Sort()
    return r

def showGraphCollection(grColl,outName,xtitle,ytitle,yran=None):
    """shows a collection of graphs in a canvas with a tidy format"""

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.06)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.02)
    
    gr=ROOT.TMultiGraph()
    leg=ROOT.TLegend(0.15,0.93,0.95,0.88)
    leg.SetNColumns(len(grColl))

    for sample,grList in grColl:

        for ig in xrange(0,len(grList)):
            g=grList[ig]
            gr.Add(g,'lX')
            if ig>0: continue
            leg.AddEntry(g,g.GetTitle(),'l')

    gr.Draw('a')
    gr.GetXaxis().SetTitleSize(0.05)
    gr.GetYaxis().SetTitleSize(0.05)
    gr.GetYaxis().SetTitleOffset(1.1)
    gr.GetXaxis().SetTitle(xtitle)
    gr.GetYaxis().SetTitle(ytitle)
    gr.GetYaxis().SetRangeUser(gr.GetYaxis().GetXmin(),1.2*gr.GetYaxis().GetXmax())
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.Draw()
    if yran:
        gr.GetYaxis().SetRangeUser(yran[0],yran[1])

    c.SetGridx()
    c.SetGridy()

    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.15,0.96,'#bf{CMS} #it{simulation preliminary}')

    c.Modified()
    c.Update()    
    for ext in ['png','pdf']:
        c.SaveAs('{0}.{1}'.format(outName,ext))
    c.Delete()

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)
for flav in ['b','c','udsg']:
    for p in ['','jrest_']:
        allGr=[]
        allRatioGr=[]
        for sample,tlist,color in [('nominal',[''],'#000000'),
                                   ('CR',['erdON','gluonMove','qcdBased'],'#bf5b17'),
                                   ('FSR',['fsrdn','fsrup'],'#beaed4'),
                                   ('ISR',['isrdn','isrup'],'#386cb0'),
                                   ('UE',['uedn','ueup'],'#7fc97f'),
                                   ('hdamp',['hdampdn','hdampup'],'#f0027f')]: 
            sgr=[]
            srgr=[]
            color=ROOT.TColor.GetColor(color)
            for t in tlist:
                pfix='' if len(t)==0 else '_'+t
                url='${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/era2016/expTageff%s.root'%pfix
                url=ROOT.gSystem.ExpandPathName(url)
                fIn=ROOT.TFile.Open(url)
                sgr.append(fIn.Get('%s%s'%(p,flav)))                
                sgr[-1].SetTitle(sample)
                sgr[-1].SetMarkerColor(color)
                sgr[-1].SetLineColor(color)
                sgr[-1].SetLineWidth(2)
                try:
                    srgr.append( getRatio( sgr[-1], allGr[0][1][0],'%s2nominal'%t ) )
                except Exception as e:
                    print e
                    pass
                fIn.Close()
            allGr.append( (sample,sgr) )
            allRatioGr.append( (sample,srgr) )

        yran=None if len(p)==0 else (0.95,1.05)
        ytitle='Tagging efficiency' if len(p)==0 else 'p_{T,rec}/p_{T,gen}'
        showGraphCollection(allGr,'%s%s'%(p,flav),xtitle='Jet p_{T} [GeV]',ytitle=ytitle,yran=yran)
        showGraphCollection(allRatioGr,'%s%s_ratio'%(p,flav),xtitle='Jet p_{T} [GeV]',ytitle=ytitle+' ratio to nominal',yran=yran)
