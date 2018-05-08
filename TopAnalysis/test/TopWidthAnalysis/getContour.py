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
parser.add_option("--wids", type="string", dest="wids"  , default="20,40,50,60,70,80,90,100,110,120,130,140,150,160,180,200,220,240,260,280,300,350,400")
parser.add_option("--extraText", type="string", dest="exText" , default="",     help="additional text to plot")

(options, args) = parser.parse_args()


wids   = options.wids.split(',')
masses = options.mass.split(',')
TwoDim = options.mass != "" and len(masses) > 1

c    = ROOT.TCanvas("","",600,400)

# create base arrays for eventual tgraph
nPoints=len(wids)
x    = ROOT.TVector(nPoints)
y    = ROOT.TVector(nPoints)
z    = ROOT.TVector(nPoints)

# Utility Function: Make 2D graph from array of triplets
def makeGraph2D(x) :
    xVec,yVec,zVec = array('d'),array('d'),array('d')
    for iG in range(0,len(x)) :
        xVec.append(x[iG][0])
        yVec.append(x[iG][1])
        zVec.append(x[iG][2])
    gr=ROOT.TGraph2D(len(x),xVec,yVec,zVec)
    gr.SetTitle("")
    return gr

# two dimensional scan if >1 mass
def get2DContour() :

    fOut = ROOT.TFile("2DInfo_%s.root"%options.outnm,"RECREATE")

    grPoints=[]
    for i,j in [(a,b) for a in xrange(0,len(wids)) for b in xrange(0,len(masses))]:
        if masses[j]=="1728" and wids[i]=="100" : continue

        fIn0 = ROOT.TFile(options.indir+"/hypotest_100vs"+wids[i]+'_m1725vs'+masses[j]+options.extname+"/testStat_scan0n_PL.root")
        fIn1 = ROOT.TFile(options.indir+"/hypotest_100vs"+wids[i]+'_m1725vs'+masses[j]+options.extname+"/testStat_scan1n_PL.root")
        tree0 = fIn0.Get("limit")
        tree1 = fIn1.Get("limit")

        if not tree0 or not tree1 :
            print "Skipping missing file",wids[i],masses[j]
        else :
            tree0.Draw("limit","quantileExpected==-1")
            tree1.Draw("limit","quantileExpected==-1")

            x     = 1.324*float(wids[i])/100
            y     = float(masses[j])/10
            z     = ( -tree0.GetV1()[0] + tree1.GetV1()[0])
            grPoints+=[(x,y,z)]

        fIn0.Close()
        fIn1.Close()

    # create graphs
    minWid=min([a[0] for a in grPoints])
    maxWid=max([a[0] for a in grPoints])
    minMass=min([a[1] for a in grPoints])
    maxMass=max([a[1] for a in grPoints])
    h=ROOT.TH2D("h","h",300,minWid,maxWid,100,minMass,maxMass)
    contours = array('d',[0.01,2.30,4.61,5.99,6.18,9.21,11.83])
    h.SetDirectory(0)

    pullGraph=makeGraph2D(grPoints)
    pullGraph.SetTitle("")
    pullGraph.SetName("chiSqGr")
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
    xmin=ROOT.TMath.MinElement(contGrs[1].GetN(),contGrs[1].GetX())
    xmax=ROOT.TMath.MaxElement(contGrs[1].GetN(),contGrs[1].GetX())
    ymin=ROOT.TMath.MinElement(contGrs[1].GetN(),contGrs[1].GetY())
    ymax=ROOT.TMath.MaxElement(contGrs[1].GetN(),contGrs[1].GetY())
    h.GetXaxis().SetRangeUser(0.5*xmin,1.4*xmax)

    centralX =ROOT.TMath.MinElement(contGrs[0].GetN(),contGrs[0].GetX())
    centralX+=ROOT.TMath.MaxElement(contGrs[0].GetN(),contGrs[0].GetX())
    centralX/=2
    centralY =ROOT.TMath.MinElement(contGrs[0].GetN(),contGrs[0].GetY())
    centralY+=ROOT.TMath.MaxElement(contGrs[0].GetN(),contGrs[0].GetY())
    centralY/=2

    print "mtop:",centralY,"+",ymax-centralY,"/ -",centralY-ymin
    print "gamm:",centralX,"+",xmax-centralX,"/ -",centralX-xmin

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
    ROOT.gStyle.SetOptStat(0)

    ltx = ROOT.TLatex(0.6,0.935,options.exText)
    ltx.SetNDC(ROOT.kTRUE)
    ltx.SetTextSize(0.05)
    ltx.Draw()

    # save plots
    c.SaveAs(options.outdir+"%s.pdf"%options.outnm)
    c.SaveAs(options.outdir+"%s.png"%options.outnm)

    # Perform a 2D fit to a paraboloid
    paraboloid = ROOT.TF2("paraboloid","pow(((x-[3])*TMath::Cos([0]) + (y-[4])*TMath::Sin([0]))/[1],2) + pow(((x-[3])*TMath::Cos([0]) - (y-[4])*TMath::Sin([0]))/[2],2)",xmin-0.2,xmax+0.2,ymin-0.2,ymax+0.2)
    paraboloid.SetParameters(0.02,centralX-xmin,centralY-ymin,0,0);
    paraboloid.FixParameter(3,centralX)
    paraboloid.FixParameter(4,centralY)
    paraboloid.SetContour(len(contours),contours)
    h.Fit(paraboloid,"","")

    h.Draw("CONT1")
    c.SaveAs(options.outdir+"%s_fit.pdf"%options.outnm)
    c.SaveAs(options.outdir+"%s_fit.png"%options.outnm)

    fOut.cd()
    pullGraph.Write()
    h.Write()
    paraboloid.Write()
    fOut.Close()


# one dimensional contour if only one mass
def get1DContour() :
    # initialize standard arrays
    for i in xrange(0,nPoints) :
        fIn0 = ROOT.TFile(options.indir+"/hypotest_100vs"+wids[i]+options.extname+"/testStat_scan0n_PL.root")
        fIn1 = ROOT.TFile(options.indir+"/hypotest_100vs"+wids[i]+options.extname+"/testStat_scan1n_PL.root")
        tree0 = fIn0.Get("limit")
        tree1 = fIn1.Get("limit")

        #print wids[i]
        tree0.Draw("limit","quantileExpected==-1")
        tree1.Draw("limit","quantileExpected==-1")

        x[i]     = 1.324*float(wids[i])/100
        y[i]     = -tree0.GetV1()[0] + tree1.GetV1()[0]

        #if wids[i]=="100" :
        #    y[i]*=-1

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

    widGuess = 1.33
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
    minNLL=pullGraph.Eval(widGuess,0,"S")

    upLim = widGuess
    foundUp = False
    dnLim = widGuess
    foundDn = False
    print "Computing upper limit"
    while not foundUp :
        if pullGraph.Eval(upLim,0,"S") >= minNLL+1.0: foundUp = True
        else : upLim += widStep
        if upLim > 20 :
            print "Couldn't find upper limit"
            break
    print "Computing lower limit"
    while not foundDn :
        if pullGraph.Eval(dnLim,0,"S") >= minNLL+1.0: foundDn = True
        else : dnLim -= widStep
        if dnLim < -5 :
            print "Couldn't find lower limit"
            break

    print "Width: ",widGuess,"  | (",dnLim,",",upLim,")"
    print widGuess,",",dnLim,",",upLim

    # remake graph
    ynew = ROOT.TVector(nPoints)
    for i in xrange(0,nPoints) :
        ynew[i] = y[i] - minNLL
    pullGr = ROOT.TGraph(x,ynew);
    pullGr.SetFillColor(ROOT.kBlack)
    pullGr.SetLineColor(ROOT.kBlack)
    pullGr.SetMarkerStyle(20)
    pullGr.SetMarkerSize(1)
    pullGr.SetTitle("")
    pullGr.Draw("ALP")
    print minNLL,pullGr.Eval(widGuess,0,"S")

    # set the bin and axis labels
    xax=pullGr.GetXaxis()
    xax.SetTitle("Top quark decay width, #Gamma_{t} (GeV)")
    xax.SetTitleOffset(0.9);
    xax.SetRangeUser(0,2*upLim)
    yax=pullGr.GetYaxis()
    yax.SetTitle("-2 ln(L_{alt.} / L_{SM})")
    yax.SetTitleOffset(1.0);
    yax.SetLabelSize(0.035);
    yax.SetRangeUser(minNLL,minNLL+10)

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
    pullGr.SaveAs(options.outdir+"%s_plot.root"%options.outnm)

if TwoDim :
    get2DContour()
else :
    get1DContour()
