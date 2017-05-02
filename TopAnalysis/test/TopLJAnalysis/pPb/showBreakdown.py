import ROOT
import sys

"""
display the results of the fit
"""
def showBreakdown(c,total,cor,wro):

    if total.Integral()==0 : return

    stackH=ROOT.THStack('stack','')
    stackH.Add(wro, 'hist' )
    stackH.Add(cor, 'hist' )

    frame=total.Clone('frame')
    frame.Reset('ICE')
    frame.GetYaxis().SetTitleOffset(1.4)
    frame.GetYaxis().SetRangeUser(0,total.GetMaximum()*1.5)
    frame.Draw()
    stackH.Draw('histsame')
    total.Draw('same')

    leg = ROOT.TLegend(0.14,0.88,0.5,0.75)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.AddEntry(total,'total','l')
    leg.AddEntry(cor,'correct','f')
    leg.AddEntry(wro,'wrong','f')
    leg.Draw()


    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.04)
    label.DrawLatex(0.15,0.9,'#bf{CMS} #it{preliminary}')

    p=total.GetName()
    tag='#mu+jets'
    if '1f'   in p : tag ='non-iso '+tag
    if '4j'   in p : tag='#geq4j,'
    if '2q'   in p : tag+='=0b'
    if '1b1q' in p : tag+='=1b'
    if '2b'   in p : tag+='#geq2b'
    label.DrawLatex(0.68,0.9,'#scale[0.8]{#it{%s}}'%tag)

    label.DrawLatex(0.68,0.85,'#scale[0.8]{Purity: %3.2f}'%(cor.Integral(0,cor.GetNbinsX()+1)/total.Integral(0,total.GetNbinsX()+1)))


    c.cd()
    c.Modified()
    c.Update()
    #c.SaveAs('%s/%s.png'%(outDir,p))
    raw_input()

    stackH.Delete()
    frame.Delete()


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.04)
c.SetBottomMargin(0.12)

fIn=ROOT.TFile.Open(sys.argv[1])
for cat in ['1l4j2q','1l4j1b1q','1l4j2b']:
    for dist in ['m_wjj','m_bwjj','m_bwl']:
        total=fIn.Get('%s_%s'%(dist,cat))
        cor=fIn.Get('%s_%s_cor'%(dist,cat))
        ci=ROOT.TColor.GetColor('#4575b4')
        cor.SetFillColor(ci)
        cor.SetFillStyle(1001)
        wro=fIn.Get('%s_%s_wro'%(dist,cat))
        ci=ROOT.TColor.GetColor('#91bfdb')
        wro.SetFillColor(ci)
        wro.SetFillStyle(1001)
        showBreakdown(c,total,cor,wro)
