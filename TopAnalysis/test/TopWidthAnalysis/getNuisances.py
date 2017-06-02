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
parser.add_option("-i",    type="string", dest="indir"  , default="./",         help="combine output file")
parser.add_option("-n",    type="string", dest="outnm"  , default="nuisances",  help="the filename for the plot")
parser.add_option("-o",    type="string", dest="outdir" , default="./",         help="the base filename for the plot")
parser.add_option("--extraText", type="string", dest="exText" , default="",     help="additional text to plot")

(options, args) = parser.parse_args()

fIn  = ROOT.TFile(options.indir)
ws   = fIn.Get("w")
snap = ws.getSnapshot("HybridNew_mc_s__snapshot")
c    = ROOT.TCanvas("","",1000,400)

# create base arrays for eventual tgraph
nPoints = snap.getSize()

x    = ROOT.TVector(nPoints)
ex   = ROOT.TVector(nPoints)
exs  = ROOT.TVector(nPoints)
y    = ROOT.TVector(nPoints)
z    = ROOT.TVector(nPoints)
ey   = ROOT.TVector(nPoints) # nuisance pulls
sg1  = ROOT.TVector(nPoints) # background graphs
sg2  = ROOT.TVector(nPoints)

xAxTitles=[]

# initialize standard arrays
for i in xrange(0,nPoints) :
    x[i]     = 0.25 + 0.5*i
    ex[i]    = 0.3
    exs[i]   = 0.25
    y[i]     = 0
    z[i]     = 0
    sg1[i]   = 1
    sg2[i]   = 2

# get pull information and errors
it  = snap.createIterator()
var = it.Next()
while var :
    xAxTitles+=[var.GetName()]
    var = it.Next()

for i in xrange(0,len(xAxTitles)) :
    y[i]  = snap.find(xAxTitles[i]).getValV()
    ey[i] = snap.find(xAxTitles[i]).getError()

# create graphs
pullGraph = ROOT.TGraphErrors(x,y,exs,ey);
sig1Graph = ROOT.TGraphErrors(x,z,ex,sg1);
sig2Graph = ROOT.TGraphErrors(x,z,ex,sg2);

# create canvas
c.cd()
setTDRStyle()
c.SetRightMargin(0.03)
c.SetLeftMargin(0.06)
c.SetBottomMargin(0.25)
c.SetGrid(0,1)
c.cd()

# format all graphs: color
pullGraph.SetFillColor(ROOT.kBlack)
pullGraph.SetLineColor(ROOT.kBlack)
pullGraph.SetMarkerStyle(20)
pullGraph.SetMarkerSize(1)

sig1Graph.SetFillColor(15)
sig2Graph.SetFillColor(17)
sig1Graph.SetLineColor(ROOT.kNone)
sig2Graph.SetLineColor(ROOT.kNone)

# draw as a multigraph
totalGraph=ROOT.TMultiGraph()
totalGraph.Add(sig2Graph)
totalGraph.Add(sig1Graph)
totalGraph.Add(pullGraph,'p')
totalGraph.Draw("a2")

# set the bin and axis labels
xax=totalGraph.GetXaxis()
xax.SetTitle("")
i=0
for label in xAxTitles :
    bin_index = xax.FindBin(0.25+0.5*i)
    xax.SetBinLabel(bin_index,label)
    i+=1
xax.SetRangeUser(0.088,0.5+0.5*(nPoints-1)-0.088)
xax.SetNdivisions(0)

yax=totalGraph.GetYaxis()
yax.SetTitle("Pull (#sigma)")
yax.SetTitleOffset(0.4);

CMS_lumi.relPosX = 0.090
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumiTextSize = 0.55
CMS_lumi.extraOverCmsTextSize=0.90
CMS_lumi.CMS_lumi(c,4,0)

ltx = ROOT.TLatex(0.6,0.935,options.exText)
ltx.SetNDC(ROOT.kTRUE)
ltx.SetTextSize(0.05)
ltx.Draw()

# save plots
c.Modified()
c.Update()
c.SaveAs(options.outdir+"%s.pdf"%options.outnm)
c.SaveAs(options.outdir+"%s.png"%options.outnm)
