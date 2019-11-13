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

hchi2=ROOT.TH1F('hchi2',';#chi^{2};Toys',100,0,1.5*chi2_obs)
htoy=hdata.Clone('htoy')
for t in range(1,1000):
    htoy.Reset('ICE')
    nevts=ROOT.gRandom.Poisson(hdata.Integral())
    for i in range(nevts):
        htoy.Fill(hbkg.GetRandom())
    hchi2.Fill( htoy.Chi2Test(hbkg,'UWCHI2') )

hchi2.Draw('hist')
hchi2.GetYaxis().SetRangeUser(0,hchi2.GetMaximum()*1.2)
l=ROOT.TLine()
l.SetLineWidth(2)
l.SetLineColor(ROOT.kBlue)
l.DrawLine(chi2_obs,0,chi2_obs,hchi2.GetMaximum())

raw_input()
