#!/usr/bin/env python
import re
from sys import argv
from array import array
import os.path
import ROOT

def interact():
    import code
    code.InteractiveConsole(locals=globals()).interact()

ROOT.gROOT.SetBatch(True)

from pprint import pprint
from optparse import OptionParser
from tdrStyle import setTDRStyle
import CMS_lumi
parser = OptionParser()
parser.add_option("-i",    type="string", dest="indir"  , default="/afs/cern.ch/work/e/ecoleman/TOP-17-010-final/datacards_inc_scan/",         help="combine output file")
parser.add_option("-e",    type="string", dest="extname", default="_100pseudodata",         help="combine output file")
parser.add_option("-n",    type="string", dest="outnm"  , default="contournosyst",    help="the filename for the plot")
parser.add_option("-o",    type="string", dest="outdir" , default="./",         help="the base filename for the plot")
parser.add_option("--mass", type="string", dest="mass"  , default="")
parser.add_option("--wids", type="string", dest="wids"  , default="20,40,50,60,70,80,90,100,110,120,130,160,200,220,240,280,300,350,400")
parser.add_option("--extraText", type="string", dest="exText" , default="",     help="additional text to plot")

(options, args) = parser.parse_args()


wids = options.wids.split(',')
mass = options.mass.split(',')
massMap={"172.5": "",
        "":      "",
        "166.5": "simmeq166p5",
        "169.5": "simmeq169p5",
        "171.5": "simmeq171p5",
        "173.5": "simmeq173p5",
        "175.5": "simmeq175p5",
        "178.5": "simmeq178p5" }
TwoDim= options.mass != "" and len(mass) > 1

c    = ROOT.TCanvas("","",600,400)

# create base arrays for eventual tgraph
nPoints=len(wids)
x    = ROOT.TVector(nPoints)
y    = ROOT.TVector(nPoints)
z    = ROOT.TVector(nPoints)

# two dimensional scan if >1 mass
def get2DContour() :
    # create graphs
    nPoints = len(wids) * len(mass)
    pullGraph = ROOT.TGraph2D(nPoints);

    h=ROOT.TH2D("","",300,0.25*1.324,3.8*1.324,100,172.44,172.54)
    contours = array('d',[2.30,4.61,5.99,6.18,9.21,11.83])

    h.SetDirectory(0)

    iPoint=0
    for i in xrange(0,len(wids)) :
        for tmas in mass :
            fIn0 = ROOT.TFile(options.indir+"/hypotest_100vs"+wids[i]+massMap[tmas]+options.extname+"/testStat_scan0n_TEV.root")
            fIn1 = ROOT.TFile(options.indir+"/hypotest_100vs"+wids[i]+massMap[tmas]+options.extname+"/testStat_scan1n_TEV.root")
            if fIn1 is None or fIn0 is None :
                iPoint += 1
                continue

            tree0 = fIn0.Get("limit")
            tree1 = fIn1.Get("limit")

            tree0.Draw("limit","quantileExpected==-1")
            tree1.Draw("limit","quantileExpected==-1")

            x     = 1.324*float(wids[i])/100
            y     = float(tmas)
            z     = ( tree0.GetV1()[0] - tree1.GetV1()[0] )
            pullGraph.SetPoint(iPoint,x,y,z)

            #print x,y,z
            fIn0.Close()
            fIn1.Close()
            iPoint += 1

    pullGraph.SetTitle("")
    pullGraph.SetHistogram(h)
    h.SetTitle("")

    # create canvas
    c.cd()
    setTDRStyle()
    c.SetGrid(1,1)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.14)
    c.SetBottomMargin(0.125)
    c.cd()

    # fill histo
    pullGraph.Draw("SURF")
    ROOT.gPad.Update()
    h.Draw("COLZ")

    # set minimum to 0
    minVal=h.GetMinimum()
    for xbin in xrange(1,h.GetNbinsX()+1):
        for ybin in xrange(1,h.GetNbinsY()+1):
            h.SetBinContent(xbin,ybin,h.GetBinContent(xbin,ybin)-minVal)
    h.GetZaxis().SetRangeUser(0,20)

    # set contours
    h.SetContour(len(contours),contours)
    h.Draw("CONTZ LIST")
    ROOT.gPad.Update()

    contList = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    contGrs  = []
    h.Draw("COLZ")
    for i in range(0,len(contours)) :
        x = contList.At(i).At(0).Clone()
        contGrs.append(x)
        x.SetLineWidth(2)
        x.Draw('same l')
    xmin=ROOT.TMath.MinElement(contGrs[0].GetN(),contGrs[0].GetX())
    xmax=ROOT.TMath.MaxElement(contGrs[0].GetN(),contGrs[0].GetX())
    ymin=ROOT.TMath.MinElement(contGrs[0].GetN(),contGrs[0].GetY())
    ymax=ROOT.TMath.MaxElement(contGrs[0].GetN(),contGrs[0].GetY())
    h.GetXaxis().SetRangeUser(0.6*xmin,1.4*xmax)

    xx, yy, zz = ROOT.Long(0), ROOT.Long(0), ROOT.Long(0)
    h.GetMinimumBin(xx,yy,zz)
    print h.GetXaxis().GetBinCenter(xx),h.GetYaxis().GetBinCenter(yy)

    h.GetXaxis().SetTitle("Top quark width, #Gamma_{t} [GeV/c^{2}]")
    h.GetYaxis().SetTitle("Top quark mass, m_{t} [GeV/c^{2}]")
    h.GetZaxis().SetTitle("-2 #Delta log(L_{alt.}/L_{ref.})")
    h.GetXaxis().SetTitleOffset(1.15)
    h.GetYaxis().SetTitleOffset(1.275)
    h.GetZaxis().SetTitleOffset(0.9)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetZaxis().SetTitleSize(0.05)

    #interact()

    # format and print
    CMS_lumi.relPosX = 0.17
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


# one dimensional contour if only one mass
def get1DContour() :
    # initialize standard arrays
    for i in xrange(0,nPoints) :
        fIn0 = ROOT.TFile(options.indir+"/hypotest_100vs"+wids[i]+massMap[options.mass]+options.extname+"/testStat_scan0n_PL.root")
        fIn1 = ROOT.TFile(options.indir+"/hypotest_100vs"+wids[i]+massMap[options.mass]+options.extname+"/testStat_scan1n_PL.root")
        tree0 = fIn0.Get("limit")
        tree1 = fIn1.Get("limit")

        tree0.Draw("limit","quantileExpected==-1")
        tree1.Draw("limit","quantileExpected==-1")

        x[i]     = 1.324*float(wids[i])/100
        y[i]     = -tree0.GetV1()[0] + tree1.GetV1()[0]
        #print x[i],y[i]
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
    pullGraph.Draw("ALP")

    widGuess = 1.324
    widStep  = 0.0001
    foundMin = False
    print "Computing minimum"
    while not foundMin :
        nomVal = pullGraph.Eval(widGuess,0,"S")
        upVal  = pullGraph.Eval(widGuess+widStep,0,"S")
        dnVal  = pullGraph.Eval(widGuess-widStep,0,"S")

        if upVal >= nomVal and dnVal >= nomVal : foundMin = True

        if upVal <= nomVal : widGuess += widStep
        elif dnVal <= nomVal : widGuess -= widStep

    upLim = 1.324
    foundUp = False
    dnLim = 1.324
    foundDn = False
    print "Computing upper limit"
    while not foundUp :
        if pullGraph.Eval(upLim,0,"S") >= 1.0: foundUp = True
        else : upLim += widStep
    print "Computing lower limit"
    while not foundDn :
        if pullGraph.Eval(dnLim,0,"S") >= 1.0: foundDn = True
        else : dnLim -= widStep
        if dnLim < -5 :
            print "Couldn't find lower limit"
            break

    print "Width: ",widGuess,"  | (",dnLim,",",upLim,")"
    print widGuess,",",dnLim,",",upLim


    # set the bin and axis labels
    xax=pullGraph.GetXaxis()
    xax.SetTitle("Top quark decay width, #Gamma_{t} (GeV)")
    xax.SetTitleOffset(0.9);
    xax.SetRangeUser(0,2*upLim)
    yax=pullGraph.GetYaxis()
    yax.SetTitle("-2 ln(L_{alt.} / L_{SM})")
    yax.SetTitleOffset(1.0);
    yax.SetLabelSize(0.035);
    yax.SetRangeUser(0,10)

    # format and print
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
    c.SaveAs(options.outdir+"%s.root"%options.outnm)
    pullGraph.SaveAs(options.outdir+"%s_plot.root"%options.outnm)

if TwoDim :
    get2DContour(0,0,0,0)
else :
    get1DContour()
