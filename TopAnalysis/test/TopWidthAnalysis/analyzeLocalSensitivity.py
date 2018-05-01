#!/usr/bin/env python

import ROOT
import pickle

def getLocalSensitivity(hlist):
    ws=ROOT.RooWorkspace('sigws')
    ws.factory('x[0,200]')

    bin1=ROOT.RooBinning(len(hlist)-1,0,len(hlist)-1)
    bin2=ROOT.RooBinning(len(hlist[0])-1,0,len(hlist[0])-1)
    refGrid=ROOT.RooMomentMorphND.Grid(bin1,bin2)

    for i in xrange(0,len(hlist)):
        for j in xrange(0,len(hlist[i])):
            name=hlist[i][j].GetName()+'_%d_%d'%(i,j)
            data=ROOT.RooDataHist(name,name,ROOT.RooArgList(ws.var("x")),hlist[i][j])
            pdf=ROOT.RooHistPdf(name+"_pdf",name+"_pdf",ROOT.RooArgSet(ws.var("x")),data)
            getattr(ws,'import')(pdf,ROOT.RooCmdArg())

            #add pdf to the grid
            print 'Adding',pdf.GetName(),'@ (',i,j,')'
            refGrid.addPdf(ws.pdf(pdf.GetName()),i,j)

    ws.factory('alpha[0,5]')
    ws.factory('beta[0,5]')
    pdf=ROOT.RooMomentMorphND('widmorphpdf','widmorphpdf',
                              ROOT.RooArgList( ws.var('alpha'), ws.var('beta') ),
                              ROOT.RooArgList( ws.var('x') ),
                              refGrid,
                              ROOT.RooMomentMorphND.Linear)
    pdf.useHorizontalMorphing(False)
    getattr(ws,'import')(pdf,ROOT.RooCmdArg())

    #central value
    ws.var('alpha').setVal(2.0)
    ws.var('beta').setVal(2.0)

    dpdfdalpha=pdf.derivative(ws.var('alpha'))
    dpdfdbeta=pdf.derivative(ws.var('beta'))

    lsSensitivities={'alpha':ROOT.TGraph(), 'beta':ROOT.TGraph()}
    xmin,xmax=hlist[0][0].GetXaxis().GetXmin(),hlist[0][0].GetXaxis().GetXmax()
    nx=hlist[0][0].GetXaxis().GetNbins()
    dx=(xmax-xmin)/nx
    norm={'alpha':0.,'beta':0}
    for i in xrange(0,nx):
        x=xmin+i*dx
        ws.var('x').setVal(x)
        f=pdf.getVal()
        if f==0 : continue
        dfda=dpdfdalpha.getVal()
        sa=(1/f)*(dfda**2)
        norm['alpha']+=sa*dx
        dfdb=dpdfdbeta.getVal()
        sb=(1/f)*(dfdb**2)
        norm['beta']+=sb*dx
        np=lsSensitivities['alpha'].GetN()
        lsSensitivities['alpha'].SetPoint(np,x,sa)
        lsSensitivities['beta'].SetPoint(np,x,sb)
    x,y=ROOT.Double(0),ROOT.Double(0)
    print norm
    for par in lsSensitivities:        
        for i in xrange(0,lsSensitivities[par].GetN()):
            lsSensitivities[par].GetPoint(i,x,y)
            lsSensitivities[par].SetPoint(i,x,y/norm[par])
        lsSensitivities[par].Sort()
    return lsSensitivities


ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

cat='EM2blowpt'
#cat='EM2bhighpt'
#cat='EM1blowpt'
#cat='EM1bhighpt'

mwidGrid=[]

for f in ['analysis/MC13TeV_TTJets_m166v5.root',
          'analysis/MC13TeV_TTJets_m169v5.root',
          'analysis/MC13TeV_TTJets.root',
          'analysis/MC13TeV_TTJets_m175v5.root',
          'analysis/MC13TeV_TTJets_m178v5.root']:
    fIn=ROOT.TFile.Open(f)
    mwidGrid.append( [fIn.Get('%s_incmlb_w%s'%(cat,w)).Clone() for w in ['60','80','100','120','140']] )
    for h in mwidGrid[-1] :
        h.Scale(1./h.Integral())
        h.SetDirectory(0)
    fIn.Close()
print  mwidGrid
mwidSens=getLocalSensitivity(mwidGrid)


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c= ROOT.TCanvas("c","c",500,500)
c.SetTopMargin(0.06)
c.SetRightMargin(0.03)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.1)

h=mwidGrid[2][2]
h.Draw('hist')
h.GetYaxis().SetTitle('PDF')
h.GetXaxis().SetTitle('Lepton-b invariant mass [GeV]')
h.GetYaxis().SetTitleSize(0.04)
h.GetXaxis().SetTitleSize(0.04)
h.GetYaxis().SetLabelSize(0.04)
h.GetXaxis().SetLabelSize(0.04)
h.GetYaxis().SetTitleOffset(1.4)
h.SetLineWidth(2)
h.SetLineColor(1)
h.GetYaxis().SetRangeUser(0,h.GetMaximum()*1.5)

mwidSens['alpha'].SetLineColor(ROOT.kRed+1)
mwidSens['alpha'].SetFillStyle(0)
mwidSens['alpha'].SetLineWidth(2)
mwidSens['alpha'].SetLineStyle(9)
mwidSens['alpha'].Draw('l')
mwidSens['beta'].SetLineColor(ROOT.kAzure-3)
mwidSens['beta'].SetLineStyle(9)
mwidSens['beta'].SetFillStyle(0)
mwidSens['beta'].SetLineWidth(2)
mwidSens['beta'].Draw('l')

leg=ROOT.TLegend(0.15,0.93,0.6,0.69)
leg.SetTextFont(42)
leg.SetTextSize(0.03)
leg.SetFillStyle(0)
leg.SetFillColor(0)
leg.SetBorderSize(0)
leg.AddEntry(h,'PDF ( m_{t} = 172.5 GeV, #Gamma_{t} = 1.34 GeV )','l')
leg.AddEntry(mwidSens['alpha'],'#frac{1}{PDF} [ #frac{dPDF}{d(#Deltam_{t} = 3.0 GeV)} ]^{2}','l')
leg.AddEntry(mwidSens['beta'],'#frac{1}{PDF} [ #frac{dPDF}{d(#Delta#Gamma_{t} = 0.27 GeV)} ]^{2}','l')
leg.Draw()

tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
tex.DrawLatex(0.85,0.96,'13 TeV')

c.Modified()
c.Update()
for ext in ['png','pdf']:
    c.SaveAs('localsens_%s.%s'%(cat,ext))
