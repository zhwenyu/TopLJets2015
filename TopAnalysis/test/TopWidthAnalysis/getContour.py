#!/usr/bin/env python
import re
from sys import argv
import os.path
import ROOT

ROOT.gROOT.SetBatch(True)

from pprint import pprint
from optparse import OptionParser
from tdrStyle import setTDRStyle
import CMS_lumi
parser = OptionParser()
parser.add_option("-i",    type="string", dest="indir"  , default="/afs/cern.ch/work/e/ecoleman/TOP-17-010/datacards_inc_scan/",         help="combine output file")
parser.add_option("-e",    type="string", dest="extname", default="_100pseudodata",         help="combine output file")
parser.add_option("-n",    type="string", dest="outnm"  , default="contournosyst",    help="the filename for the plot")
parser.add_option("-o",    type="string", dest="outdir" , default="./",         help="the base filename for the plot")
parser.add_option("--wids", type="string", dest="wids"  , default="20,40,60,80,100,120,140,160,200,220,240,280,300,350,400")
parser.add_option("--extraText", type="string", dest="exText" , default="",     help="additional text to plot")

(options, args) = parser.parse_args()

wids = options.wids.split(',')
c    = ROOT.TCanvas("","",600,400)

# create base arrays for eventual tgraph
nPoints=len(wids)
x    = ROOT.TVector(nPoints)
y    = ROOT.TVector(nPoints)

# initialize standard arrays
for i in xrange(0,nPoints) :
    fIn0 = ROOT.TFile(options.indir+"/hypotest_100vs"+wids[i]+options.extname+"/higgsCombinescan0n.HybridNew.mH172.5.123456.root")
    fIn1 = ROOT.TFile(options.indir+"/hypotest_100vs"+wids[i]+options.extname+"/higgsCombinescan1n.HybridNew.mH172.5.123456.root")
    tree0 = fIn0.Get("limit")
    tree1 = fIn1.Get("limit")

    tree0.Draw("limit","quantileExpected==-1")
    tree1.Draw("limit","quantileExpected==-1")

    x[i]     = 1.324*float(wids[i])/100
    y[i]     = tree0.GetV1()[0] - tree1.GetV1()[0]
    print x[i],y[i]
    fIn0.Close()
    fIn1.Close()

# create graphs
pullGraph = ROOT.TGraph(x,y);

# create canvas
c.cd()
setTDRStyle()
c.SetGrid(1,1)
c.SetRightMargin(0.01)
c.SetLeftMargin(0.125)
c.SetBottomMargin(0.125)
c.cd()

# format all graphs: color
pullGraph.SetFillColor(ROOT.kBlack)
pullGraph.SetLineColor(ROOT.kBlack)
pullGraph.SetMarkerStyle(20)
pullGraph.SetMarkerSize(1)
pullGraph.SetTitle("")
pullGraph.Draw("AL")

# set the bin and axis labels
xax=pullGraph.GetXaxis()
xax.SetTitle("Top width (GeV)")
xax.SetTitleOffset(0.9);

yax=pullGraph.GetYaxis()
yax.SetTitle("-ln(L_{alt.} / L_{SM})")
yax.SetTitleOffset(1.0);
yax.SetLabelSize(0.035);

CMS_lumi.relPosX = 0.165
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumiTextSize = 0.55
CMS_lumi.extraOverCmsTextSize=0.80
CMS_lumi.CMS_lumi(c,4,0)

ltx = ROOT.TLatex(0.6,0.935,options.exText)
ltx.SetNDC(ROOT.kTRUE)
ltx.SetTextSize(0.05)
ltx.Draw()

# save plots
c.SaveAs(options.outdir+"%s.pdf"%options.outnm)
c.SaveAs(options.outdir+"%s.png"%options.outnm)
