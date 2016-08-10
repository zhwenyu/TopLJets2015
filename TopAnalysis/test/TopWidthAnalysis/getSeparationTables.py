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
parser = OptionParser(
    usage="%prog [options]",
    epilog="Collects quantiles information from signal statistics output and turns it into a nice TGraph and LaTeX table. Format of .txt files is stats__<wid>_<dist>.txt"
    )
parser.add_option("-i",    type="string", dest="indir"  , default="./",     help="directory to look for stats files in")
parser.add_option("-o",    type="string", dest="outdir" , default="./",     help="the base filename for the quantiles plot")
parser.add_option("--dist",type="string", dest="dist"   , default="mlb",    help="the observable distribution to look at")
parser.add_option("--wid", type="string", dest="widList", default="0p5w,1p0w,1p5w,2p0w,2p5w,3p0w,3p5w,4p0w,4p5w,5p0w",
        help="a list of widths to look for in stats filenames")
parser.add_option("--recreateOut", dest="recreate", default=False, action="store_true")
parser.add_option("--doAll", dest="doAll", default=False, action="store_true")

(options, args) = parser.parse_args()

# get lists to loop over
rawWidList=options.widList.split(',')
rawWidList.remove('1p0w')

rawSepGraphY = []

distToTitle={
        "mlb": ("Inclusive M_{lb}",0.73,0.83),
        "minmlb": ("Minimum M_{lb}",0.725,0.83),
        "mdrmlb": ("#DeltaR-Filtered M_{lb}",0.7,0.83),
        "incmlb": ("Inclusive M_{lb}",0.73,0.83),
        "sncmlb": ("Semi-Inclusive M_{lb}",0.65,0.83),
        "mt2mlb": ("M_{T2}^{lb} Strategy",0.73,0.83)
       }

# initialize stats
statList={
        "Separation": {},
        "$P(q_{\\rm null}>q_{\\rm alt}^{\\rm median})$": {},
        "$P(q_{\\rm alt}<q_{\\rm null}^{\\rm median})$": {},
        "CL$_s$^{\\rm exp.}$": {},
        "CL$_s$^{\\rm obs.}$": {},
        }

canvas=ROOT.TCanvas("","",800,600)
canvas.cd()

if not options.doAll :
    for wid,stat in [(w,s) for w in rawWidList for s in statList] :
        statList[stat][wid.replace('p','.').replace('w','$\\times\\Gamma_{\\rm SM}$')] = "-1"

    # loop over widths, parse array info
    for wid in rawWidList :
        latexWid = wid.replace('p','.').replace('w','$\\times\\Gamma_{\\rm SM}$')

        statsFileName="%s/stats__%s__%s.txt"%(options.indir,wid,options.dist)
        for line in open(statsFileName,"r"):
            if "separation" in line :
               statList["Separation"][latexWid]=line.split('#')[0]
            if "sep" in line and not "separation" in line:
                rawSepGraphY += [float(line.split('#')[0])]
            elif "null exceeded density" in line :
               statList["$P(q_{\\rm null}>q_{\\rm alt}^{\\rm median})$"][latexWid]=line.split('#')[0]
            elif "alt exceeded density" in line :
               statList["$P(q_{\\rm alt}<q_{\\rm null}^{\\rm median})$"][latexWid]=line.split('#')[0]
            elif "cls expected" in line :
               statList["CL$_s$^{\\rm exp.}$"][latexWid]=line.split('#')[0]
            elif "cls observed" in line :
               statList["CL$_s$^{\\rm obs.}$"][latexWid]=line.split('#')[0]
            else : continue

    tabletex=open("%s/separationTable__%s.tex"%(options.outdir,options.dist), 'w')

    tabletex.write("\\begin{table}\n")
    tabletex.write("\\centering\n")
    tabletex.write("\\begin{tabular}{|l%s|} \\hline "%(len(rawWidList)*"|c"))

    for wid in rawWidList :
        last = (wid == rawWidList[-1])
        toWrite = wid.replace('p','.').replace('w','$\\times\\Gamma_{\\rm SM}$')
        tabletex.write("%12s %s"%(toWrite, " \\\\\\hline\n" if last else " & "))

    # oh python, why do you not maintain order in maps?
    newStatList = sorted([key for key in statList])[::-1]
    for stat in newStatList :
        tabletex.write("%30s "%(stat))
        for wid in sorted(statList[stat]) :
            value = statList[stat][wid]
            if value == "-1" : value = ""
            tabletex.write("& %12s "%(value))

        if not stat == newStatList[-1] :
            tabletex.write(" \\\\\\hline")
        tabletex.write("\n")
    tabletex.write("\\end{tabular}\n")
    tabletex.write("\\caption{}\n")
    tabletex.write("\\label{tab:separation%s}\n"%options.dist)
    tabletex.write("\\end{table}\n")
    tabletex.close()

    setTDRStyle()
    xarr=[1.324*float(wid[0:3].replace('p','.').replace('w','')) for wid in sorted([key for key in statList["CL$_s$^{\\rm exp.}$"]])]
    yarr=[float(statList["CL$_s$^{\\rm exp.}$"][wid]) for wid in sorted([key for key in statList["CL$_s$^{\\rm exp.}$"]])]
    x=ROOT.TVector(len(rawWidList))
    y=ROOT.TVector(len(rawWidList))

    for ix in xrange(0,len(xarr)):
        x[ix] = float(xarr[ix])

    for iy in xrange(0,len(yarr)):
        y[iy] = float(yarr[iy])

    clsGr=ROOT.TGraph(x,y)
    clsGr.GetXaxis().SetTitle("Generator-level #Gamma_{t} [GeV]")
    clsGr.GetYaxis().SetTitle("CL_{s}^{exp.}")
    clsGr.SetTitle("")
    clsGr.SetName("CLsExp_%s"%options.dist)

    sepGraphY=ROOT.TVector(len(yarr))
    for iy in xrange(0,len(yarr)) :
        sepGraphY[iy] = rawSepGraphY[iy]

    sepGr=ROOT.TGraph(x,sepGraphY)
    sepGr.GetXaxis().SetTitle("Generator-level #Gamma_{t} [GeV]")
    sepGr.GetYaxis().SetTitle("Separation")
    sepGr.SetTitle("")
    sepGr.SetName("Separation_%s"%options.dist)


    #tpad1=ROOT.TPad("tpad1","",0,0,1,1)
    #tpad2=ROOT.TPad("tpad2","",0,0,1,1)
    #tpad2.SetFillStyle(4000)

    ymin=0
    ymax=2000
    dy=(ymax-ymin)/0.8
    xmin=-3
    xmax=3
    dx=(xmax-xmin)/0.8

    clsGr.Draw("ALP")

    fout=ROOT.TFile("%s/statsPlots.root"%(options.outdir),"UPDATE" if not options.recreate else "RECREATE")
    fout.cd()
    sepGr.Write()
    clsGr.Write()
    fout.Close()

    # draw dist information if available
    if options.dist in [key for key in distToTitle] :
        DistInfo,xpos,ypos=distToTitle[options.dist]
        DistLaTeX=ROOT.TLatex(xpos,ypos, DistInfo)
        DistLaTeX.SetNDC(ROOT.kTRUE)
        DistLaTeX.SetTextSize(0.04)
        DistLaTeX.Draw()

        TMassInfo=ROOT.TLatex(.745,.77,"m_{t} = 172.5 GeV")
        TMassInfo.SetNDC(ROOT.kTRUE)
        TMassInfo.SetTextSize(0.03)
        TMassInfo.Draw()

    CMS_lumi.relPosX = 0.180
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.extraOverCmsTextSize=0.50
    CMS_lumi.lumiTextSize = 0.55
    CMS_lumi.CMS_lumi(canvas,4,0)

    canvas.SetGrid(True)
    #canvas.SetLogy(True)
    canvas.SaveAs("%s/CLsPlot__%s.png"%(options.outdir,options.dist))
    canvas.SaveAs("%s/CLsPlot__%s.pdf"%(options.outdir,options.dist))

if options.doAll :

    graphs={}

    for dist in options.dist.split(',') :

        for wid,stat in [(w,s) for w in rawWidList for s in statList] :
            statList[stat][wid.replace('p','.').replace('w','$\\times\\Gamma_{\\rm SM}$')] = "-1"

        # loop over widths, parse array info
        for wid in rawWidList :
            latexWid = wid.replace('p','.').replace('w','$\\times\\Gamma_{\\rm SM}$')

            statsFileName="%s/stats__%s__%s.txt"%(options.indir,wid,dist)
            for line in open(statsFileName,"r"):
                if "separation" in line :
                   statList["Separation"][latexWid]=line.split('#')[0]
                if "sep" in line and not "separation" in line:
                    rawSepGraphY += [float(line.split('#')[0])]
                elif "null exceeded density" in line :
                   statList["$P(q_{\\rm null}>q_{\\rm alt}^{\\rm median})$"][latexWid]=line.split('#')[0]
                elif "alt exceeded density" in line :
                   statList["$P(q_{\\rm alt}<q_{\\rm null}^{\\rm median})$"][latexWid]=line.split('#')[0]
                elif "cls expected" in line :
                   statList["CL$_s$^{\\rm exp.}$"][latexWid]=line.split('#')[0]
                elif "cls observed" in line :
                   statList["CL$_s$^{\\rm obs.}$"][latexWid]=line.split('#')[0]
                else : continue

        tabletex=open("%s/separationTable__%s.tex"%(options.outdir,dist), 'w')

        tabletex.write("\\begin{table}\n")
        tabletex.write("\\centering\n")
        tabletex.write("\\begin{tabular}{|l%s|} \\hline "%(len(rawWidList)*"|c"))

        for wid in rawWidList :
            last = (wid == rawWidList[-1])
            toWrite = wid.replace('p','.').replace('w','$\\times\\Gamma_{\\rm SM}$')
            tabletex.write("%12s %s"%(toWrite, " \\\\\\hline\n" if last else " & "))

        # oh python, why do you not maintain order in maps?
        newStatList = sorted([key for key in statList])[::-1]
        for stat in newStatList :
            tabletex.write("%30s "%(stat))
            for wid in sorted(statList[stat]) :
                value = statList[stat][wid]
                if value == "-1" : value = ""
                tabletex.write("& %12s "%(value))

            if not stat == newStatList[-1] :
                tabletex.write(" \\\\\\hline")
            tabletex.write("\n")
        tabletex.write("\\end{tabular}\n")
        tabletex.write("\\caption{}\n")
        tabletex.write("\\label{tab:separation%s}\n"%dist)
        tabletex.write("\\end{table}\n")
        tabletex.close()

        xarr=[1.324*float(wid[0:3].replace('p','.').replace('w','')) for wid in sorted([key for key in statList["CL$_s$^{\\rm exp.}$"]])]
        yarr=[float(statList["CL$_s$^{\\rm exp.}$"][wid]) for wid in sorted([key for key in statList["CL$_s$^{\\rm exp.}$"]])]
        x=ROOT.TVector(len(rawWidList))
        y=ROOT.TVector(len(rawWidList))

        for ix in xrange(0,len(xarr)):
            x[ix] = float(xarr[ix])

        for iy in xrange(0,len(yarr)):
            y[iy] = float(yarr[iy])

        clsGr=ROOT.TGraph(x,y)
        clsGr.GetXaxis().SetTitle("Generator-level #Gamma_{t} [GeV]")
        clsGr.GetYaxis().SetTitle("CL_{s}^{exp.}")
        clsGr.SetTitle("")
        clsGr.SetName("CLsExp_%s"%dist)

        sepGraphY=ROOT.TVector(len(yarr))
        for iy in xrange(0,len(yarr)) :
            sepGraphY[iy] = rawSepGraphY[iy]

        sepGr=ROOT.TGraph(x,sepGraphY)
        sepGr.GetXaxis().SetTitle("Generator-level #Gamma_{t} [GeV]")
        sepGr.GetYaxis().SetTitle("Separation")
        sepGr.SetTitle("")
        sepGr.SetName("Separation_%s"%dist)

        graphs["Separation_%s"%dist]=sepGr
        graphs["CLsExp_%s"%dist]=clsGr

    setTDRStyle()

    kOrangeList=[0,1,2,4,3]

    i=0
    clsLeg=ROOT.TLegend(.60,.85,.85,.60)
    for dist in options.dist.split(',') :
        tgr=graphs["CLsExp_%s"%dist]
        tgr.SetLineColor(ROOT.kOrange+kOrangeList[i])
        tgr.SetMarkerColor(ROOT.kOrange+kOrangeList[i])
        tgr.SetMarkerStyle(ROOT.kFullCircle)
        clsLeg.AddEntry(tgr,distToTitle[dist][0],"LP")

        tf=ROOT.TF1("ClsFit_%s"%dist,"expo",1.324,1.324*3.5)
        tf.SetLineColor(ROOT.kOrange+kOrangeList[i])
        tgr.Fit(tf,"NQ","",1.324,1.324*3.5)


        if i==0:
            tgr.GetXaxis().SetRangeUser(1,4)
            tgr.GetYaxis().SetRangeUser(0,1)
        tgr.Draw("|P" if i>0 else "AP")
        tf.Draw("SAME")

        i+=1
    clsLeg.Draw()

    ROOT.gStyle.SetOptFit(0)

    l=ROOT.TLine()
    l.SetLineColor(ROOT.kRed)
    l.DrawLine(1,0.05,4.08,0.05)

    CMS_lumi.relPosX=0.18
    CMS_lumi.extraText="Simulation Preliminary"
    CMS_lumi.extraOverCmsTextSize=.5
    CMS_lumi.lumiTextSize=.7
    CMS_lumi.CMS_lumi(canvas,4,0)
    canvas.SaveAs("%s/CLsPlot.png"%(options.outdir))
    canvas.SaveAs("%s/CLsPlot.pdf"%(options.outdir))


    canvas.Clear()

    i=0
    sepLeg=ROOT.TLegend(.65,.15,.85,.45)
    for dist in options.dist.split(',') :
        tgr=graphs["Separation_%s"%dist]
        tgr.SetLineColor(ROOT.kRed+i)
        tgr.SetMarkerStyle(ROOT.kNone)

        clsLeg.AddEntry(tgr,distToTitle[dist][0],"L")

        tgr.Draw("SAME AL" if i>0 else "AL")

        i+=1


    CMS_lumi.CMS_lumi(canvas,4,0)
    canvas.SaveAs("%s/SeparationPlot.png"%(options.outdir))
    canvas.SaveAs("%s/SeparationPlot.pdf"%(options.outdir))

