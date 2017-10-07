import ROOT
import sys

COLORS=[ROOT.kGray+3,ROOT.kBlue-3, ROOT.kRed-4, ROOT.kOrange-3]
WIDTH=[4,3,3,1]

"""
"""
def doPlot(fileList,plot='mjj'):

    histos={}
    for i in xrange(0,len(fileList)):
        tag,url=fileList[i].split(':')
        fIn=ROOT.TFile.Open(url)

        for cat in ['1l4j2q','1l4j1b1q','1l4j2b']:
            for t in ['','_cor','_wro']:
                
                key=(cat,t)
                if not key in histos : histos[key]=[]
                try:
                    histos[key].append( fIn.Get('%s_e%s%s'%(plot,cat,t)).Clone('%s_%s%s_%d'%(plot,cat,t,i)) )
                    histos[key][-1].Add( fIn.Get('%s_mu%s%s'%(plot,cat,t) ) )
                    histos[key][-1].SetDirectory(0)
                    histos[key][-1].SetTitle(tag)
                    histos[key][-1].SetLineWidth(WIDTH[i-1])
                    histos[key][-1].SetLineColor(COLORS[i-1])
                    histos[key][-1].SetMarkerColor(COLORS[i-1])
                    histos[key][-1].SetMarkerStyle(1)
                    if '#DeltaR' in tag:
                        histos[key][-1].SetFillColor(ROOT.kOrange-2)
                        histos[key][-1].SetFillStyle(1001)
                    if 'leading' in tag:
                        histos[key][-1].SetFillColor(ROOT.kGray+2)
                        histos[key][-1].SetFillStyle(3354)
                except:
                    pass
                    
                    
        fIn.Close()


    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)

    for key in histos:
        if len(histos[key])==0 : continue
        drawOpt='hist'
        for h in histos[key]:
            h.Scale(1./h.Integral())
            h.Draw(drawOpt)
            h.GetYaxis().SetRangeUser(0,h.GetMaximum()*2.3)
            h.GetYaxis().SetTitleOffset(1.2)
            drawOpt='histsame'
            h.GetYaxis().SetTitleOffset(1.2)
            h.GetXaxis().SetTitleOffset(1.0)
            h.GetYaxis().SetTitleSize(0.045)
            h.GetXaxis().SetTitleSize(0.045)
            h.GetYaxis().SetLabelSize(0.035)
            h.GetXaxis().SetLabelSize(0.035)
            title=h.GetXaxis().GetTitle()
            h.GetYaxis().SetTitle('Events (a.u.)')
            h.SetLineWidth(3)
            if 'W_{jj}' in title : h.GetXaxis().SetTitle("m_{jj'} [GeV]")

        leg=c.BuildLegend(0.7,0.95,0.95,0.7)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)

        label = ROOT.TLatex()
        label.SetNDC()
        label.SetTextFont(42)
        label.SetTextSize(0.04)
        label.DrawLatex(0.15,0.9,'#scale[1.2]{#bf{CMS}} #scale[0.75]{#it{Simulation Supplementary}}')
        label.DrawLatex(0.645,0.965,'#scale[0.9]{pPb (#sqrt{s_{NN}} = 8.16 TeV)}')
        tag='e^{#pm} / #mu^{#pm} + #geq4j'
        if '1f'   in key[0] : tag ='non-iso '+tag
        if '2q'   in key[0] : tag+=' (=0b)'
        if '1b1q' in key[0] : tag+=' (=1b)' 
        if '2b'   in key[0] : tag+=' (#geq2b)'
        label.DrawLatex(0.15,0.85,'#scale[0.9]{#bf{%s}}'%tag)
        print tag
        assign='inc. '
        if 'cor' in key[1]: assign='t#bar{t} correct'
        if 'wro' in key[1]: assign='t#bar{t} wrong'
        label.DrawLatex(0.15,0.8,'#scale[0.9]{%s assignments}'%(assign))
        label.DrawLatex(0.15,0.75,'#scale[0.9]{pN#rightarrow t#bar{t}+X (PYTHIA6 Z2*)}')
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('%s_%s%s.%s'%(plot,key[0],key[1],ext))
            print '%s_%s%s.%s'%(plot,key[0],key[1],ext)

"""
"""
def main():

    for plot in ['mjj','mthad','mtlep']:
        doPlot(fileList=sys.argv[1:],plot=plot)

"""
"""
if __name__ == "__main__":
    main()
