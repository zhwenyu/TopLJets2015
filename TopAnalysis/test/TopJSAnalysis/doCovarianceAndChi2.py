import ROOT
ROOT.gROOT.SetBatch(True)
import optparse
import os,sys
import json
import re
from collections import OrderedDict
from math import sqrt
from array import *
import random
import numpy
import copy

debug = True

"""
steer the script
"""
def main():
    
    cmsLabel='#bf{CMS} #it{preliminary}'
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option(     '--mcUnc',        dest='mcUnc'  ,      help='common MC related uncertainty (e.g. lumi)',        default=0,              type=float)
    parser.add_option(     '--com',          dest='com'  ,        help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default='data/era2016/samples.json',              type='string')
    parser.add_option( '--systJson', dest='systJson', help='json with list of systematics', default='data/era2016/syst_samples.json', type='string')
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default='unfolding/result',              type='string')
    parser.add_option('', '--inDirToys',       dest='inDirToys' ,      help='input toy directory',                default='unfolding/toys',              type='string')
    parser.add_option('-O', '--outDir',      dest='outDir' ,     help='output directory',                default='unfolding/covariance',              type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=16551.,              type=float)
    parser.add_option('--obs', dest='obs',  default='mult', help='observable [default: %default]')
    parser.add_option('--flavor', dest='flavor',  default='all', help='flavor [default: %default]')
    (opt, args) = parser.parse_args()
    
    # statistical covariance
    
    toyfile = '%s/%s_charged_%s_toys.root'%(opt.inDirToys, opt.obs, opt.flavor)
    fInToy = ROOT.TFile.Open(toyfile)
    
    pseudoresults = []
    counter = 0
    while True:
        pseudoresult = []
        h = fInToy.Get('Unfolded_' + str(counter))
        if not h: break
        for i in range(1, h.GetNbinsX()+1):
            pseudoresult.append(h.GetBinContent(i)/h.GetBinWidth(i))
        integral = sum(pseudoresult)
        for i in range(len(pseudoresult)):
            pseudoresult[i] = pseudoresult[i]/integral
        pseudoresults.append(pseudoresult)
        counter += 1
    
    print('Imported', counter, 'toy experiments')
    
    x = numpy.array(pseudoresults).T
    print(x)
    statcov = numpy.cov(x)
    print(statcov)
    statcov_reduced = numpy.delete(statcov, 3, 0)
    statcov_reduced = numpy.delete(statcov_reduced, 3, 1)
    print(statcov_reduced)
    print(numpy.linalg.det(statcov_reduced))
    print(numpy.linalg.inv(statcov_reduced))
    
    #rootoutfile = ROOT.TFile.Open(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_cov.root', 'RECREATE')
    #rootoutfile.cd()
    
    # 11-class RdBu http://colorbrewer2.org/#type=diverging&scheme=RdBu&n=11
    stops = array('d', [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    red   = array('d')
    green = array('d')
    blue  = array('d')
    colors = [[103,0,31],
              [178,24,43],
              [214,96,77],
              [244,165,130],
              [253,219,199],
              [247,247,247],
              [209,229,240],
              [146,197,222],
              [67,147,195],
              [33,102,172],
              [5,48,97]]
    for color in colors:
        red.append(color[0]/255.)
        green.append(color[1]/255.)
        blue.append(color[2]/255.)
    ROOT.TColor.CreateGradientColorTable(11, stops, red[::-1], green[::-1], blue[::-1], 30)
    
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas('c','c',500,500)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    c.cd()
    
    h = fInToy.Get('Unfolded_0')
    dataStatCovNorm = ROOT.TH2D('dataStatCovNorm', '', h.GetNbinsX(), h.GetXaxis().GetXbins().GetArray(), h.GetNbinsX(), h.GetXaxis().GetXbins().GetArray())
    axistitle = h.GetXaxis().GetTitle().replace('generated ', '')
    dataStatCovNorm.SetXTitle(axistitle)
    dataStatCovNorm.SetYTitle(axistitle)
    
    for i in range(1, dataStatCovNorm.GetNbinsX()+1):
        for j in range(1, dataStatCovNorm.GetNbinsY()+1):
            dataStatCovNorm.SetBinContent(i, j, statcov[i-1][j-1])
    
    dataStatCovNorm.GetZaxis().SetRangeUser(-max(abs(statcov.min()), abs(statcov.max())),
                                             max(abs(statcov.min()), abs(statcov.max())))
    dataStatCovNorm.Draw('colz')
    
    ROOT.gPad.Update()
    tl1 = ROOT.TLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl1.Draw()
    tl2 = ROOT.TLine(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl2.Draw()
    
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.16,0.91, cmsLabel)
    txt.DrawLatex(0.63,0.97, '#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
    txt.DrawLatex(0.16,0.85, 'Statistical covariance')
    
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_cov_stat.pdf')
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_cov_stat.png')
    
    
    # Read lists of syst samples
    allSystVars = ['jec_CorrelationGroupMPFInSitu',
                   'jec_RelativeFSR',
                   'jec_CorrelationGroupUncorrelated',
                   'jec_FlavorPureGluon',
                   'jec_FlavorPureQuark',
                   'jec_FlavorPureCharm',
                   'jec_FlavorPureBottom',
                   'jer',
                   'btag_heavy',
                   'btag_light',
                   'csv_heavy',
                   'csv_light',
                   'tracking'
                  ]
    varList = []
    for var in allSystVars:
        varList.append(var+'_up')
        varList.append(var+'_down')
    varList += ['evtgen',
                'm171v5',
                'm173v5',
                'herwig',
                'isrup',
                'isrdn',
                'fsrup',
                'fsrdn',
                'hdampup',
                'hdampdn',
                'ueup',
                'uedn',
                'erdON',
                'qcdBased']

    expSystSamplesList = []
    for var in varList:
        expSystSamplesList.append(['MC13TeV_TTJets_'+var, [832., 0., '', 't#bar{t} '+var]])
    
    resultfile = '%s/%s_charged_%s_result.root'%(opt.inDir, opt.obs, opt.flavor)
    fIn=ROOT.TFile.Open(resultfile)
    
    # reference
    hnominal    = fIn.Get('MC13TeV_TTJets_Unfolded')
    nominal     = []
    for i in range(1, hnominal.GetNbinsX()+1):
        nominal.append(hnominal.GetBinContent(i))
    
    systcov = copy.copy(statcov)
    systcov -= systcov
    
    for slist in [expSystSamplesList]:
        for tag,sample in slist: 
            if tag in ['MC13TeV_TTJets_cflip']: continue
            if not 't#bar{t}' in sample[3] : continue

            hsyst = fIn.Get(tag+'_Unfolded')
            if not hsyst:
                print(tag, 'not loaded')  
                continue
            
            syst = []
            for i in range(1, hsyst.GetNbinsX()+1):
                syst.append(hsyst.GetBinContent(i))
            
            x = numpy.array([nominal, syst]).T
            cov = numpy.cov(x)
            systcov += cov
            
    dataSystCovNorm = dataStatCovNorm.Clone('dataSystCovNorm')
    dataSystCovNorm.Reset()
    
    for i in range(1, dataSystCovNorm.GetNbinsX()+1):
        for j in range(1, dataSystCovNorm.GetNbinsY()+1):
            dataSystCovNorm.SetBinContent(i, j, systcov[i-1][j-1])
    
    dataSystCovNorm.GetZaxis().SetRangeUser(-max(abs(systcov.min()), abs(systcov.max())),
                                             max(abs(systcov.min()), abs(systcov.max())))
    dataSystCovNorm.Draw('colz')
    
    ROOT.gPad.Update()
    tl1 = ROOT.TLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl1.Draw()
    tl2 = ROOT.TLine(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl2.Draw()
    
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.16,0.91, cmsLabel)
    txt.DrawLatex(0.63,0.97, '#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
    txt.DrawLatex(0.16,0.85, 'Systematic covariance')
    
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_cov_syst.pdf')
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_cov_syst.png')
    
    dataCovNorm = dataStatCovNorm.Clone('dataCovNorm')
    dataCovNorm.Reset()
    
    cov = statcov + systcov
    
    for i in range(1, dataCovNorm.GetNbinsX()+1):
        for j in range(1, dataCovNorm.GetNbinsY()+1):
            dataCovNorm.SetBinContent(i, j, cov[i-1][j-1])
    
    dataCovNorm.GetZaxis().SetRangeUser(-max(abs(cov.min()), abs(cov.max())),
                                         max(abs(cov.min()), abs(cov.max())))
    dataCovNorm.Draw('colz')
    
    ROOT.gPad.Update()
    tl1 = ROOT.TLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl1.Draw()
    tl2 = ROOT.TLine(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl2.Draw()
    
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.16,0.91, cmsLabel)
    txt.DrawLatex(0.63,0.97, '#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
    txt.DrawLatex(0.16,0.85, 'Total covariance')
    
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_cov.pdf')
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_cov.png')
    
    print(numpy.linalg.det(cov))
    cov_reduced = numpy.delete(cov, 0, 0)
    cov_reduced = numpy.delete(cov_reduced, 0, 1)
    print(cov_reduced)
    print(numpy.linalg.det(cov_reduced))
    print(numpy.linalg.inv(cov_reduced))
    
    print(opt.obs, opt.flavor, 'chi2 nominal',  returnChi2(fIn, cov_reduced, nominal, 'nominalGen'))
    print(opt.obs, opt.flavor, 'chi2 fsr up',   returnChi2(fIn, cov_reduced, nominal, 'FSRUpGen'))
    print(opt.obs, opt.flavor, 'chi2 fsr down', returnChi2(fIn, cov_reduced, nominal, 'FSRDownGen'))
    print(opt.obs, opt.flavor, 'chi2 herwig',   returnChi2(fIn, cov_reduced, nominal, 'herwigGen'))
    
def returnChi2(fIn, cov_reduced, nominal, prediction):
    hpred = fIn.Get(prediction)
    pred  = []
    for i in range(1, hpred.GetNbinsX()+1):
        pred.append(hpred.GetBinContent(i))
    diff = []
    for i in range(len(nominal)):
        diff.append(pred[i] - nominal[i])
    
    chi2 = numpy.array(diff[1:]).T.dot(numpy.linalg.inv(cov_reduced).dot(numpy.array(diff[1:])))
    ndf  = hpred.GetNbinsX()-2
    prob = ROOT.TMath.Prob(chi2, ndf)
    return chi2/ndf, prob

def normalizeAndDivideByBinWidth(hist):
    hist.Scale(1./hist.Integral())
    for i in range(1, hist.GetNbinsX()+1):
        hist.SetBinContent(i, hist.GetBinContent(i)/hist.GetBinWidth(i))
        hist.SetBinError  (i, hist.GetBinError(i)  /hist.GetBinWidth(i))
    return hist
        
"""
for execution from another script
"""
if __name__ == "__main__":
    main()
    #sys.exit(main())

