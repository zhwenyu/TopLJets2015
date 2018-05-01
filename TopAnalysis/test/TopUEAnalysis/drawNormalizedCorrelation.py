import sys
import ROOT
from UETools import getNormalizedPerColumn

pname=sys.argv[2]

fIn=ROOT.TFile.Open(sys.argv[1])
h=fIn.Get(pname)
#h.RebinX()
#h.RebinY()
hnorm=getNormalizedPerColumn(h)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
ROOT.gStyle.SetPaintTextFormat("3.0f");

c=ROOT.TCanvas('c','c',550,500)
c.SetTopMargin(0.06)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.12)
c.SetRightMargin(0.16)

hnorm.Draw('colz text')
hnorm.GetYaxis().SetTitleOffset(1.1)
hnorm.GetYaxis().SetTitleSize(0.05)
hnorm.GetXaxis().SetTitleSize(0.05)
hnorm.GetZaxis().SetTitleSize(0.05)

#px=hnorm.ProfileX()
#px.SetMarkerStyle(20)
#px.Draw('e1same')

tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.1,0.96,'#bf{CMS} #it{simulation}')
tex.DrawLatex(0.8,0.96,'#sqrt{s}=13 TeV')
print h.GetCorrelationFactor()
print hnorm.GetCorrelationFactor()
c.RedrawAxis()
c.Modified()
c.Update()
raw_input()
for ext in ['png','pdf','root']: c.SaveAs('%s.%s'%(pname,ext))
