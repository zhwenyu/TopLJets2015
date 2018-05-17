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
snap = ws.getSnapshot("MultiDimFit")
c    = ROOT.TCanvas("","",1000,400)

# create base arrays for eventual tgraph
nPoints = snap.getSize()

x    = ROOT.TVector(nPoints)
xx   = ROOT.TVector(4*nPoints)
ex   = ROOT.TVector(4*nPoints)
exs  = ROOT.TVector(nPoints)
y    = ROOT.TVector(nPoints)
z    = ROOT.TVector(4*nPoints)
ey   = ROOT.TVector(nPoints) # nuisance pulls
sg1  = ROOT.TVector(4*nPoints) # background graphs
sg2  = ROOT.TVector(4*nPoints)

xAxTitles=[]

# initialize standard arrays
for i in xrange(0,nPoints) :
    x[i]     = 0.25 + 0.5*i
    exs[i]   = 0 #0.25
    y[i]     = -1

for i in xrange(0,4*nPoints) :
    xx[i]    = 0.25 + 0.5/4 * i
    ex[i]    = 0.25 / 4 #- (0.23/100 if i % 4==0 else 0)
    z[i]     = 0
    sg1[i]   = 1
    sg2[i]   = 2


# get pull information and errors
it  = snap.createIterator()
var = it.Next()
while var :
    passes=var.GetName() != "CMS_th1x"
    passes=passes and var.GetName() != "CMS_channel"
    passes=passes and var.GetName() != "ZERO"
    passes=passes and var.GetName() != "ONE"
    passes=passes and "In" not in var.GetName()
    if passes :
        varValCrit = abs(var.getValV()) < 0.1
        varErrCrit = abs(var.getError()-1) < 0.4
        if not varValCrit or not varErrCrit :
            xAxTitles+=[var.GetName()]
    var = it.Next()

for i in xrange(0,len(xAxTitles)) :
    y[i]  = snap.find(xAxTitles[i]).getValV()
    ey[i] = snap.find(xAxTitles[i]).getError()
    if xAxTitles[i] == "mtop" :
        print "mtop:",y[i],'+/-',ey[i]

# create graphs
pullGraph = ROOT.TGraphErrors(x,y,exs,ey);
sig1Graph = ROOT.TGraphErrors(xx,z,ex,sg1);
sig2Graph = ROOT.TGraphErrors(xx,z,ex,sg2);

# create canvas
c.cd()
setTDRStyle()
c.SetRightMargin(0.01)
c.SetLeftMargin(0.06)
c.SetBottomMargin(0.3)
c.SetGrid(0,1)
c.cd()

# format all graphs: color
pullGraph.SetFillColor(ROOT.kBlack)
pullGraph.SetLineColor(ROOT.kBlack)
pullGraph.SetMarkerStyle(20)
pullGraph.SetMarkerSize(0.75)

sig1Graph.SetFillColor(42)
sig2Graph.SetFillColor(18)
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
xax.SetLabelSize(0.05)
xax.SetTitle("")
i=0
for label in xAxTitles :
    bin_index = xax.FindBin(0.25+0.5*i)
    xax.SetBinLabel(bin_index,label)
    i+=1
xax.SetRangeUser(0.088,0.25+0.5*(len(xAxTitles)-1))
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

line = ROOT.TLine(0,0,0.25+0.5*(len(xAxTitles)-1),0)
line.SetLineStyle(ROOT.kSolid)
line.SetLineColor(ROOT.kBlack)
line.Draw()

# save plots
c.Modified()
c.Update()
c.SaveAs(options.outdir+"%s.pdf"%options.outnm)
c.SaveAs(options.outdir+"%s.png"%options.outnm)
