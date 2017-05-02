#!/usr/bin/env python2.7

import ROOT
import optparse
import os,sys

from run2dFits import showFitResult,FITVARS,XRANGES

"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(False)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--output',    dest='output',    default='plots/Data8TeV_pp',          type='string',   help='output directory [%default]')
    parser.add_option('-i', '--input',     dest='input',     default='workspace_Data8TeV_pp.root', type='string',   help='workspace [%default]')
    parser.add_option('-v', '--verbose',   dest='verbose',   default=0,                            type=int,        help='Verbose mode [%default]')
    (opt, args) = parser.parse_args()

    #read workspace
    fIn=ROOT.TFile(opt.input)
    w=fIn.Get('w')

    w.loadSnapshot('fit')
    allVars=w.allVars()
    varIter = allVars.createIterator()
    var=varIter.Next()
    toFix=['sig_','eb_','w_','fqcd_','e4jtoe3j','Nbkg_1f4j1q','Nbkg_1f4j1b']
    while var :
        varName=var.GetName()
        for v in toFix:
            if varName.find(v)!=0 : continue
            w.var(varName).setConstant(True)
            if v=='eb_' or v=='e4jtoe3j' :
                w.var(varName).setVal(1.0)
                w.var(varName).setError(0.0)
            break
        var=varIter.Next()

    yieldsList=ROOT.RooArgList(w.var('Nsig'),
                        #w.var('Nbkg_1f4j1q'),
                        #w.var('Nbkg_1l4j1b'),
                        #w.var('Nbkg_1l4j1q'),
                        w.var('Nbkg_1l4j2b'))
    pdf=w.pdf('pdf')
    data=w.data('data')
    redData = data.reduce(ROOT.RooFit.Cut("sample==sample::1l4j2b"))
    pdf.fitTo(redData,ROOT.RooFit.Extended())

    for fitVar in FITVARS:
        xmin,xmax=XRANGES[fitVar]
        showFitResult(fitVar=fitVar,data=data,pdf=pdf,categs=['1l4j2b'],w=w,showComponents=['S_*','S_*,W_*'],rangeX=(xmin,xmax),postfix='2dsplot',outDir=opt.output)

    sData = ROOT.RooStats.SPlot("s"+redData.GetName(),"An SPlot",redData,pdf,yieldsList)
    sigwdata = ROOT.RooDataSet('sigwdata','sigwdata',redData,redData.get(),'','Nsig_sw')
    bkgwdata = ROOT.RooDataSet('bkgwdata','bkgwdata',redData,redData.get(),'','Nbkg_1l4j2b_sw')

    #compare background vs signal like hypothesis
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.02)
    c.SetBottomMargin(0.1)
    VARTITLES={'dr_jj':'#DeltaR(j,j)',
               'pt_wl':'p_{T}(W_{l}) [GeV]',
               'y_wl':'y(w_{lep})',
               'y_l':'y(l)',
               'pt_wjj':'p_{T}(W_{had}) [GeV]',
               'y_wjj':'y(w_{had})',
               'y_bhad':'y(b_{had})',
               'y_blep':'y(b_{lep})',
               'y_bwjj':'y(t_{had})',
               'y_bwl':'y(t_{lep})',
               'dy_bwjjbwl':'|y(t_{lep})-y(t_{had})|',
               'sy_bwjjbwl':'|y(t_{lep})+y(t_{had})|'}

    for var in VARTITLES:
        varTitle=VARTITLES[var]
        w.var(var).setBins(10)
        frame=w.var(var).frame(ROOT.RooFit.Name('sig'))
        bkgwdata.plotOn(frame,
                        ROOT.RooFit.LineColor(ROOT.kGray),
                        ROOT.RooFit.MarkerStyle(24),
                        ROOT.RooFit.MarkerColor(ROOT.kGray),
                        ROOT.RooFit.Name('bkg')
                        )
        sigwdata.plotOn(frame)
        frame.Draw()
        frame.GetXaxis().SetTitle(varTitle)
        label = ROOT.TLatex()
        label.SetNDC()
        label.SetTextFont(42)
        label.SetTextSize(0.04)
        label.DrawLatex(0.15,0.9,'#bf{CMS} #it{preliminary}')
        label.DrawLatex(0.8,0.96,'(#scale[0.8]{#sqrt{s}=8.16 TeV})')
        leg=ROOT.TLegend(0.15,0.88,0.4,0.75)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        leg.AddEntry('sig','Signal','ep')
        leg.AddEntry('bkg','Background','ep')
        leg.Draw()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('%s/%s_splot.%s'%(opt.output,var,ext))



if __name__ == "__main__":
    main()
