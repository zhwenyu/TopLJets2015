#!/usr/bin/env python
import re
from sys import argv
import os.path
import ROOT

ROOT.gROOT.SetBatch(true)

from pprint import pprint
from optparse import OptionParser
parser = OptionParser(
    usage="%prog [options]",
    epilog="Collects quantiles information from signal statistics output and turns it into a nice TGraph and LaTeX table. Format of .txt files is stats__<wid>_<dist>.txt"
    )
parser.add_option("-i",    type="string", dest="indir"  , default="./",     help="directory to look for stats files in")
parser.add_option("-o",    type="string", dest="outdir" , default="./",     help="the base filename for the quantiles plot")
parser.add_option("--dist",type="string", dest="dist"   , default="mlb",    help="the observable distribution to look at")
parser.add_option("--wid", type="string", dest="widList", default="0p5w,1p0w,1p5w,2p0w,2p5w,3p0w,3p5w,4p0w,4p5w,5p0w",
        help="a list of widths to look for in stats filenames")

(options, args) = parser.parse_args()

# get lists to loop over
rawWidList=options.widList.split(',')

# initialize stats
statList={
        "Separation": {}
        "$P(q_{\\rm null}>q_{\\rm alt}^{\\rm median})$": {},
        "$P(q_{\\rm alt}<q_{\\rm null}^{\\rm median})$": {},
        "CL$_s$^{\\rm exp.}$": {},
        "CL$_s$^{\\rm obs.}$": {},
        }

for wid,stat in [(w,s) for w in rawWidList for stat in statList] :
    statList[stat][wid.replace('p','.').replace('w','$\\times\\Gamma_{\\rm SM}$')] = "-1"

# loop over widths, parse array info
for wid in rawWidList :
    latexWid = wid.replace('p','.').replace('w','$\\times\\Gamma_{\\rm SM}$')

    statsFileName="%s/stats__%s_%s.txt"%(options.indir,wid,options.dist)
    for line in open(statsFileName,"r"):
        if "separation" in line :
           statList["Separation"][latexWid]=line.split('#')[0]
        elif "null exceeded density" in line :
           statList["$P(q_{\\rm null}>q_{\\rm alt}^{\\rm median})$"][latexWid]=line.split('#')[0]
        elif "alt exceeded density" in line :
           statList["$P(q_{\\rm alt}<q_{\\rm null}^{\\rm median})$"][latexWid]=line.split('#')[0]
        elif "cls expected" in line :
           statList["CL$_s$^{\\rm exp.}$"][latexWid]=line.split('#')[0]
        elif "cls observed" in line :
           statList["CL$_s$^{\\rm obs.}$"][latexWid]=line.split('#')[0]
        else : continue

tabletex=open("%s/separationTable__%s.tex"%(options.outdir,options.dist))

tabletex.write("\\begin{table}\n")
tabletex.write("\\centering\n")
tabletex.write("\\begin{tabular}{|l%s|} \\hline "%(len(rawWidList)*"|c"))

for wid in rawWidList
    last = (wid == rawWidList[-1])
    toWrite = wid.replace('p','.').replace('w','$\\times\\Gamma_{\\rm SM}$')
    tabletex.write("%12s %s"%(toWrite, " \\\\\\hline\n" if last else " & "

# oh python, why do you not maintain order in maps?
for stat in reversed(sorted([key for key in statList])) :
    tabletex.write("%30s "%(stat))
    for wid in sorted(statList[stat]) :
        value = statList[stat][wid]
        if value == "-1" : value = ""
        tabletex.write("& %12s "%(value))

    if not stat == reversed(sorted([key for key in statList]))[-1] :
        tabletex.write(" \\\\\\hline")
    tabletex.write("\n")
tabletex.write("\\end{tabular}\n")
tabletex.write("\\caption{}\n"%options.dist)
tabletex.write("\\label{tab:separation%s}\n"%options.dist)
tabletex.write("\\end{table}\n")
tabletex.close()

ROOT.gROOT.ProcessLine(".L tdrStyle.C")
ROOT.gROOT.ProcessLine("setTDRStyle()")

clsGr=ROOT.TGraph(len(rawWidList),
    [1.324*float(wid.replace('p','.').replace('w','')) for wid in rawWidList],
    [float(statList["CL$_s$^{\\rm exp.}$"][wid]) for wid in rawWidList])
clsGr.GetXaxis().SetTitle("Generator-level #Gamma_{t}")
clsGr.GetYaxis().SetTitle("CL_{s}^{exp.}")

canvas=ROOT.TCanvas()
canvas.cd()

# CMS text
CMSLine="CMS"
CP=ROOT.TLatex(0.12,0.92, CMSLine)
CP.SetNDC(ROOT.kTRUE)
CP.SetTextSize(0.05)
CP.Draw()

# Lumi
CMSLineLumi="#sqrt{s}=13 TeV, 2.3 fb^{-1}"
CP1=ROOT.TLatex(0.67,0.92, CMSLineLumi)
CP1.SetNDC(ROOT.kTRUE)
CP1.SetTextSize(0.04)
CP1.Draw()

# ExtraText
CMSLineExtra="#bf{#it{Preliminary}}"
CP2=ROOT.TLatex(0.195,0.92, CMSLineExtra)
CP2.SetNDC(ROOT.kTRUE)
CP2.SetTextSize(0.04)
CP2.Draw()

clsGr.Draw("AL")

canvas.SaveAs("%s/CLsPlot__%s.png"%(options.outdir,options.dist))
canvas.SaveAs("%s/CLsPlot__%s.pdf"%(options.outdir,options.dist))
