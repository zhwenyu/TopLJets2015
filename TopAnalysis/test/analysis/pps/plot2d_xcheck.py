import ROOT

def finalize(c,extra,out):
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.SetTextAlign(31)
    tex.DrawLatex(0.95,0.96,extra)
    c.Modified()
    c.Update()   
    c.SaveAs(out+'.png')


sig=ROOT.TChain('data')
for i in [120,130,140,150]:
    sig.AddFile("/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p04/Z_m_X_1200_xangle_%d_2017_preTS2_opt_v1_simu_reco.root"%i)

bkg=ROOT.TChain('data')
for i in 'BCDEF':
    bkg.AddFile("/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p04/Data13TeV_2017%s_MuonEG.root"%i)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.15)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)

sig.Draw("ypp:mpp >> h(25,500,2500,50,-1,1)","wgt*(bosonpt>40 && mixType==0 && cat==169)","goff")
h=ROOT.gDirectory.Get('h')
h.GetYaxis().SetTitle("Di-proton rapidity")
h.GetXaxis().SetTitle('Di-proton mass [GeV]')
h.Draw('colz')
finalize(c,'m_{X}=1.2 TeV (pre-TS2)','sig_yppmpp')


sig.Draw("mmiss:mpp >> h(25,500,2500,50,0,2500)","wgt*(bosonpt>40 && mixType==0 && cat==169)","goff")
h=ROOT.gDirectory.Get('h')
h.GetYaxis().SetTitle("Missing mass [GeV]")
h.GetXaxis().SetTitle('Di-proton mass [GeV]')
h.Draw('colz')
finalize(c,'m_{X}=1.2 TeV (pre-TS2)','sig_mmissmpp')

bkg.Draw("ypp:mpp >> h(25,500,2500,50,-1,1)","wgt*(mmiss>0 && bosonpt>40 && mixType==0 && cat==143)","goff")
h=ROOT.gDirectory.Get('h')
h.GetYaxis().SetTitle("Di-proton rapidity")
h.GetXaxis().SetTitle('Di-proton mass [GeV]')
h.Draw('colz')
finalize(c,'e#mu (2017)','emu_yppmpp')


bkg.Draw("mmiss:mpp >> h(25,500,2500,50,0,2500)","wgt*(mmiss>0 && bosonpt>40 && mixType==0 && cat==143)","goff")
h=ROOT.gDirectory.Get('h')
h.GetYaxis().SetTitle("Missing mass [GeV]")
h.GetXaxis().SetTitle('Di-proton mass [GeV]')
h.Draw('colz')
finalize(c,'e#mu (2017)','emu_mmissmpp')


bkg.Draw("ypp:mpp >> h(25,500,2500,50,-1,1)","wgt*(mmiss>0 && bosonpt>40 && mixType==1 && cat==143)","goff")
h=ROOT.gDirectory.Get('h')
h.GetYaxis().SetTitle("Di-proton rapidity")
h.GetXaxis().SetTitle('Di-proton mass [GeV]')
h.Draw('colz')
finalize(c,'e#mu (2017)','emubkg_yppmpp')


bkg.Draw("mmiss:mpp >> h(25,500,2500,50,0,2500)","wgt*(mmiss>0 && bosonpt>40 && mixType==1 && cat==143)","goff")
h=ROOT.gDirectory.Get('h')
h.GetYaxis().SetTitle("Missing mass [GeV]")
h.GetXaxis().SetTitle('Di-proton mass [GeV]')
h.Draw('colz')
finalize(c,'e#mu (2017)','emubkg_mmissmpp')

