import os
import ROOT

baseDir='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p04/'
t=ROOT.TChain('data')
for i in 'BCDEF':
    t.AddFile(os.path.join(baseDir,'Data13TeV_2017%s_SinglePhoton.root'%i))

t.Draw('mmiss>>hdata(50,0,2500)',    'wgt*(mixType==0 && cat==22 && xangle==150 && bosonpt>95 && mmiss>0)','goff')
hdata=ROOT.gDirectory.Get('hdata')
hdata.SetDirectory(0)
t.Draw('mmiss>>hbkg(50,0,2500)', 'wgt*(mixType==1 && cat==22 && xangle==150 && bosonpt>95 && mmiss>0)','goff')
hbkg=ROOT.gDirectory.Get('hbkg')
hbkg.SetDirectory(0)

hbkg.Scale(hdata.Integral()/hbkg.Integral())


chi2_obs=hdata.Chi2Test(hbkg,'UWCHI2')

hchi2=ROOT.TH1F('hchi2',';#chi^{2};Toys',50,20,120)
htoy=hdata.Clone('htoy')
for t in range(1,1000):
    htoy.Reset('ICE')
    nevts=ROOT.gRandom.Poisson(hdata.Integral())
    for i in range(nevts):
        htoy.Fill(hbkg.GetRandom())
    hchi2.Fill( htoy.Chi2Test(hbkg,'UWCHI2') )

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
tex.DrawLatex(0.95,0.96,'#scale[0.6]{2.3 fb^{-1} (13 TeV, 150#murad)}')
c.Modified()
c.Update()   
c.SaveAs('chi2fig25bottomleft.pdf')
c.SaveAs('chi2fig25bottomleft.png')
