import ROOT
from TopLJets2015.TopAnalysis.Plot import fixExtremities, scaleTo

def getEff(h,forward=True):

    """build an efficiency curve for a given cut"""

    eff=ROOT.TGraphErrors()
    eff.SetTitle(h.GetTitle())
    eff.SetMarkerStyle(h.GetMarkerStyle())
    eff.SetMarkerColor(h.GetMarkerColor())
    eff.SetLineColor(h.GetLineColor())
    nbins=h.GetNbinsX()
    err=ROOT.Double(0)    
    total=h.Integral()
    for xbin in xrange(1,nbins+1):        
        n=h.IntegralAndError(xbin,nbins,err)
        eff.SetPoint(xbin-1,h.GetXaxis().GetBinLowEdge(xbin),n/total)
        eff.SetPointError(xbin-1,0,float(err)/total)
    return eff

def getROC(bkgH,sigH):

    """build a ROC type of curve: eff signal vs eff background"""

    roc=ROOT.TGraphErrors()
    roc.SetTitle(sigH.GetTitle())
    roc.SetMarkerStyle(sigH.GetMarkerStyle())
    roc.SetMarkerColor(sigH.GetMarkerColor())
    roc.SetLineColor(sigH.GetLineColor())

    nbins=bkgH.GetNbinsX()
    errSig=ROOT.Double(0)
    errBkg=ROOT.Double(0)
    totalBkg=bkgH.Integral()
    totalSig=sigH.Integral()
    for xbin in xrange(1,nbins+1):
        nbkg=bkgH.IntegralAndError(xbin,nbins,errBkg)
        nsig=sigH.IntegralAndError(xbin,nbins,errSig)
        roc.SetPoint(xbin-1,nsig/totalSig,nbkg/totalBkg)
        roc.SetPointError(xbin-1,float(errSig)/totalSig,float(errBkg)/totalBkg)

    return roc



def getCMSLabel():
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')


def doEffPlot(ch,d,plotter,bkg_plotter):

    """do the efficiency plots"""
    
    #background
    bkg_gen=bkg_plotter.Get('gen{0}_{1}/gen{0}_{1}_DY'.format(ch,d) )
    bkg_rec=bkg_plotter.Get('gen{0}rec_{1}/gen{0}rec_{1}_DY'.format(ch,d) )
    bkg_rec.SetDirectory(0)
    bkg_rec.SetLineColor(1)
    bkg_rec.SetMarkerColor(1)
    bkg_rec.SetMarkerStyle(20)
    fixExtremities(bkg_gen)
    fixExtremities(bkg_rec)
    genH=[bkg_gen]
    recH=[bkg_rec]

    #signals
    gen_dir=plotter.GetDirectory('gen{0}_{1}'.format(ch,d))
    rec_dir=plotter.GetDirectory('gen{0}rec_{1}'.format(ch,d))
    colors=['#666666','#8c510a','#7fc97f', '#f0027f', '#fdc086', '#358CEF', '#bf5b17']
    for key in gen_dir.GetListOfKeys():

        ci=ROOT.TColor.GetColor(colors[len(recH)-2])
        marker=20+len(recH)/2

        genH.append(key.ReadObj())
        genH[-1].SetLineColor(ci)
        genH[-1].SetMarkerColor(ci)
        genH[-1].SetMarkerStyle(marker)
        genH[-1].SetFillStyle(0)
        genH[-1].SetLineWidth(2)
        fixExtremities(genH[-1])

        proc=key.GetName()
        proc=proc.replace('gen{0}_{1}'.format(ch,d),'')
        recH.append(rec_dir.Get('gen{0}rec_{1}{2}'.format(ch,d,proc)))
        recH[-1].SetLineColor(ci)
        recH[-1].SetMarkerColor(ci)
        recH[-1].SetMarkerStyle(marker)
        recH[-1].SetFillStyle(0)
        recH[-1].SetLineWidth(2)
        fixExtremities(recH[-1])
    
    #compute efficiency
    effGr=[]
    accEffGr=[]
    rocGr=[]
    for i in xrange(0,len(genH)):
        accEffGr.append( ROOT.TGraphAsymmErrors() )
        accEffGr[-1].SetName('acceff%d'%i)
        accEffGr[-1].SetTitle(recH[i].GetTitle())
        accEffGr[-1].SetMarkerStyle(recH[i].GetMarkerStyle())
        accEffGr[-1].SetMarkerColor(recH[i].GetMarkerColor())
        accEffGr[-1].SetLineColor(recH[i].GetLineColor())
        accEffGr[-1].BayesDivide(recH[i],genH[i])

        effGr.append( getEff(recH[i] ) )
        effGr[-1].SetName('eff%d'%i)

        if i==0: continue
        rocGr.append( getROC( recH[0],recH[i] ) )
        rocGr[-1].SetName('roc%d'%i)

    
    #show plots
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.02)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetGridx()
    c.SetGridy()

    for item in [(accEffGr,recH[0].GetXaxis().GetTitle(),            'Acc x #varepsilon',            False, False, None, (0,1),    'tr', '%s_%s_acceff'%(ch,d)),
                 (effGr,   'Cut on ' + recH[0].GetXaxis().GetTitle(),'#varepsilon',                  False, False, None, (0,1),    'br', '%s_%s_eff'%(ch,d)),
                 (rocGr,   'Signal Acc x #varepsilon',               'Background Acc x #varepsilon', False, True,  None, (1e-3,1), 'tl', '%s_%s_roc'%(ch,d))
                 ]:

        grColl,xtit,ytit,logx,logy,xran,yran,legPos,outName=item

        c.SetLogx(logx)
        c.SetLogy(logy)
        if legPos=='tl':
            leg=ROOT.TLegend(0.15,0.95,0.45,0.95-0.04*len(grColl))
        elif legPos=='tr':
            leg=ROOT.TLegend(0.65,0.95,0.95,0.95-0.04*len(grColl))
        else:
            leg=ROOT.TLegend(0.65,0.45,0.95,0.45-0.04*len(grColl))
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        mg=ROOT.TMultiGraph()
        for g in grColl:
            mg.Add(g,'p')
            leg.AddEntry(g,g.GetTitle(),'ep')
        mg.Draw('a')
        if xtit:
            mg.GetXaxis().SetTitle(xtit)
        if xran:
            mg.GetYaxis().SetRangeUser(*xran)
        if ytit:
            mg.GetYaxis().SetTitle(ytit)
        if yran:
            mg.GetYaxis().SetRangeUser(*yran)
        leg.Draw()
        getCMSLabel()
        c.Modified()
        c.Update()   
        c.RedrawAxis() 
        c.SaveAs(outName+'.png')
    
    c.Clear()
    c.SetLogy(True)
    c.SetLogx(False)
    leg=ROOT.TLegend(0.65,0.95,0.95,0.95-0.04*len(genH))
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    drawOpt='hist'
    for h in genH:
        scaleTo(h,1.0)
        h.Draw(drawOpt)
        h.GetYaxis().SetRangeUser(1e-4,1.0)
        h.GetYaxis().SetTitle('PDF')
        leg.AddEntry(h,h.GetTitle(),'lf')
        drawOpt='histsame'
    leg.Draw()
    getCMSLabel()
    c.Modified()
    c.Update()
    c.RedrawAxis()
    c.SaveAs('%s_%s.png'%(ch,d))


bkg_plotter=ROOT.TFile.Open('plots/zx_sel/bkg_plotter.root')
plotter=ROOT.TFile.Open('plots/zx_sel/plotter.root')
dists=['ptll','mll','drll']
chs=['ee','mm','em']
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)
for ch in chs:
    for d in dists:
        doEffPlot(ch,d,plotter,bkg_plotter)

