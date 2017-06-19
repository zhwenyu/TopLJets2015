#!/usr/bin/env python

import ROOT
import sys
import os
import optparse

"""
"""
def saveCanvas(c,name,plotExt=['png','pdf']):
    for ext in plotExt:
        c.SaveAs('%s.%s'%(name,ext))

"""
"""
def summarizeResultsFromToys(url,outDir='./'):

    tag=os.path.splitext(os.path.basename(url))[0]
    fIn=ROOT.TFile.Open(url)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    for cname in ['cdatasub','cdata','cmig','cunfolded','cdatafold']:
        c=fIn.Get(cname)
        c.Draw()
        saveCanvas(c=c,name='%s/%s_%s'%(outDir,cname,tag))

    cscan=fIn.Get('cscan')
    cscan.Draw()
    saveCanvas(c=cscan,name='%s/tauscan_%s'%(outDir,tag))

    for hname in ['global_bias','global_pulls']:
        h=fIn.Get(hname)
        cscan.Clear()
        h.Draw('e1')
        h.SetMarkerStyle(20)
        h.GetYaxis().SetTitle('Toys x bins')
        h.Fit('gaus','MLQ+')
        f=h.GetFunction('gaus')

        txt=ROOT.TLatex()
        txt.SetNDC()
        txt.SetTextFont(42)
        txt.SetTextSize(0.04)
        txt.DrawLatex(0.15,0.90,'#bf{CMS} simulation preliminary')
        txt.DrawLatex(0.85,0.96,'#scale[0.8]{(13 TeV)}')
        txt.DrawLatex(0.15,0.85,'#mu=%3.2f#pm%3.2f'%(f.GetParameter(1),f.GetParError(1)))
        txt.DrawLatex(0.15,0.80,'#sigma=%3.2f#pm%3.2f'%(f.GetParameter(2),f.GetParError(2)))
        txt.DrawLatex(0.15,0.75,'Mean=%3.2f#pm%3.2f'%(h.GetMean(),h.GetMeanError()))
        txt.DrawLatex(0.15,0.70,'RMS=%3.2f#pm%3.2f'%(h.GetRMS(),h.GetRMSError()))

        cscan.Modified()
        cscan.Update()
        saveCanvas(c=cscan,name='%s/%s_%s'%(outDir,tag,hname))


    cscan.SetRightMargin(0.15)
    cscan.SetLeftMargin(0.15)
    for hname in ['bin_bias','bin_pulls']:
        h=fIn.Get(hname)
        cscan.Clear()

        gr=ROOT.TGraphErrors()
        gr.SetMarkerStyle(20)
        for xbin in xrange(1,h.GetNbinsX()+1):
            hproj=h.ProjectionY("px",xbin,xbin)
            hproj.Fit('gaus','MLQ+')
            f=hproj.GetFunction('gaus')
            mu=f.GetParameter(1)
            sigma=f.GetParameter(2)
            np=gr.GetN()
            gr.SetPoint(np,h.GetXaxis().GetBinCenter(xbin),mu)
            gr.SetPointError(np,0,sigma)
            hproj.Delete()

        h.Draw('colz')
        gr.Draw('p')
        h.RebinY()
        h.GetYaxis().SetTitleOffset(1.5)
        h.GetZaxis().SetTitle('Toys')
        txt=ROOT.TLatex()
        txt.SetNDC()
        txt.SetTextFont(42)
        txt.SetTextSize(0.04)
        txt.DrawLatex(0.2,0.90,'#bf{CMS} simulation preliminary')
        txt.DrawLatex(0.75,0.96,'#scale[0.8]{(13 TeV)}')
        ROOT.gStyle.SetPalette(ROOT.kAquamarine)
        
        cscan.Modified()
        cscan.Update()
        saveCanvas(c=cscan,name='%s/%s_%s'%(outDir,tag,hname))

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',     dest='input',      help='input file [%default]',  default=None,                    type='string')
    (opt, args) = parser.parse_args()

    summarizeResultsFromToys(url=opt.input,outDir=os.path.dirname(opt.input))

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

