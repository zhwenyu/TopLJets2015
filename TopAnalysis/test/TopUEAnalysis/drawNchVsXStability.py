import ROOT

#plot='nchvsnvtx'
plot='nchvsrho'
fIn=ROOT.TFile.Open('UEanalysis/plots/plotter.root')
data=fIn.Get('%s_EM/%s_EM'%(plot,plot))
ttbar=fIn.Get('%s_EM/%s_EM_t#bar{t}'%(plot,plot))

profData=data.ProfileY("dataprof")
profData.SetMarkerStyle(20)
profttbar=ttbar.ProfileY("ttbarprof")
profttbar.SetMarkerStyle(24)
profttbar.SetMarkerColor(1)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetBottomMargin(0.1)
c.SetLeftMargin(0.12)
c.SetTopMargin(0.05)
c.SetRightMargin(0.02)

profttbar.Draw('e1')
profData.Draw('e1same')
if 'vsnvtx' in plot : 
    profttbar.GetXaxis().SetRangeUser(10,40)
    profttbar.Fit('pol1','RQ+','same',10,40)
    profData.Fit('pol1','RQ+','same',10,40)
if 'vsrho' in plot : 
    profttbar.GetXaxis().SetRangeUser(10,35)
    profttbar.Fit('pol1','RQ+','same',10,35)
    profData.Fit('pol1','RQ+','same',10,35)

profttbar.GetYaxis().SetRangeUser(18,28)
profttbar.GetXaxis().SetTitle('Vertex multiplicity' if 'vsnvtx' in plot else '#rho')
profttbar.GetYaxis().SetTitle('Charged particle multiplicity')

txt=ROOT.TLatex()
txt.SetNDC(True)
txt.SetTextFont(42)
txt.SetTextSize(0.1)
txt.SetTextAlign(12)
txt.SetTextSize(0.04)
txt.DrawLatex(0.15,0.9,'#bf{CMS} #it{preliminary}')
txt.DrawLatex(0.15,0.85,'#it{reconstruction level}')
txt.DrawLatex(0.75,0.97,'#scale[0.8]{36.5 fb^{-1} (13 TeV)}')

leg = ROOT.TLegend(0.55,0.4,0.95,0.3)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
varT='N_{vtx}' if 'vsnvtx' in plot else '#rho'
leg.AddEntry(profData,'data N_{ch}/%s=%3.2f#pm%3.2f'%(varT,profData.GetFunction('pol1').GetParameter(1),profData.GetFunction('pol1').GetParError(1)),'pe')
leg.AddEntry(profttbar,'t#bar{t} N_{ch}/%s=%3.2f#pm%3.2f'%(varT,profttbar.GetFunction('pol1').GetParameter(1),profttbar.GetFunction('pol1').GetParError(1)),'p')
leg.Draw()

c.Modified()
c.Update()
for ext in ['png','pdf']:
    c.SaveAs('%s.%s'%(plot,ext))

