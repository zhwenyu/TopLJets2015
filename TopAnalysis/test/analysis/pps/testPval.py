import os
import ROOT

url='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/analysis_1exc/bkg_ptll0/plotter_embkg.root'
d='mmiss_0' #multi-multi
d='mmiss_6' #low pu
d='mmiss_15' #nch<15 

#get bkg and data from file
fIn=ROOT.TFile.Open(url)
hbkg=fIn.Get('{}/mmiss_bkg_MuonEG_obs_background'.format(d))
hbkg.SetDirectory(0)
hdata=fIn.Get('{}/mmiss_data_MuonEG_obs'.format(d))
hdata.SetDirectory(0)
fIn.Close()

#scale to te same yields
hbkg.Scale(hdata.Integral()/hbkg.Integral())

#observed chi^2
chi2_obs=hdata.Chi2Test(hbkg,'UWCHI2')

#throw toys and build a distribution of expected chi^2
hchi2=ROOT.TH1F('hchi2',';#chi^{2};Toys',50,20,120)
htoy=hdata.Clone('htoy')
for t in range(1,2000):
    htoy.Reset('ICE')
    nevts=ROOT.gRandom.Poisson(hdata.Integral())
    for i in range(nevts):
        htoy.Fill(hbkg.GetRandom())
    hchi2.Fill( htoy.Chi2Test(hbkg,'UWCHI2') )


import numpy as np

gr=ROOT.TGraph()
for i,p in enumerate([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]):
    prob=np.array([p])
    q=np.zeros_like(prob)
    hchi2.GetQuantiles(1,q,prob)
    gr.SetPoint(i,q[0],p)

pval=gr.Eval(chi2_obs)
if pval>0.5: pval=1-pval
gr.Print('all')


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)

hchi2.Draw('hist')
hchi2.GetYaxis().SetRangeUser(0,hchi2.GetMaximum()*1.2)
l=ROOT.TLine()
l.SetLineWidth(2)
l.SetLineColor(ROOT.kBlue)
l.DrawLine(chi2_obs,0,chi2_obs,hchi2.GetMaximum())

tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
tex.SetTextAlign(31)
tex.DrawLatex(0.95,0.96,'#scale[0.7]{e#mu, 37.1 fb^{-1} (13 TeV)}')

tex.DrawLatex(0.9,0.9,'p-val:{:3.3f}'.format(pval))
c.Modified()
c.Update()   
c.SaveAs('chi2_{}.pdf'.format(d))
c.SaveAs('chi2_{}.png'.format(d))


