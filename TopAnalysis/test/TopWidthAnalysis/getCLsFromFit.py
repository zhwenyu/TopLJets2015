#!/usr/bin/env python
import re
from sys import argv
import os.path
import ROOT

#ROOT.gROOT.SetBatch(True)

from pprint import pprint
from optparse import OptionParser
from tdrStyle import setTDRStyle
import CMS_lumi
parser = OptionParser(
    usage="%prog [options]",
    epilog="Collects confidence level information from CLs fits"
    )
parser.add_option("-i",     type="string", dest="indir"  , default="./",     help="directory to look for stats files in")
parser.add_option("--iname",type="string", dest="iname"  , default="statsPlots.root", help="filename containing fit info")
parser.add_option("--dist", type="string", dest="dist"   , default="incmlb",    help="the observable distribution to look at")
parser.add_option("--nomWid", type="string", dest="inwid"   , default="1.324",    help="the observable distribution to look at")
parser.add_option("--prep", type="string", dest="prepost"   , default="post",    help="the observable distribution to look at")
parser.add_option("--lims", type="string", dest="limList"   , default="0.05,0.01",    help="the observable distribution to look at")
parser.add_option("--unblind", dest="unblind", default=False, action='store_true',  help="the observable distribution to look at")
parser.add_option("--splineVersion", dest="splineV", type="string", default='2',  help="the tspline type to use")

(options, args) = parser.parse_args()


def getSplineIntersection(yvalue, spline, scanRes=0.01, startValue=float(options.inwid),minValue=0,maxValue=6,tolerance=0.0005) :

    scanX=startValue
    lastPoint=float('inf')
    upperLimit=None
    lowerLimit=None

    currentScanRes=scanRes
    scanDirection=1
    startAbove=spline.Eval(startValue) >= yvalue

    # Get the upper limit
    while (scanX < maxValue) :
        splineVal=spline.Eval(scanX)

        # check if we've crossed the intersection yvalue
        switchCondition=(startAbove != (splineVal >= yvalue))

        if switchCondition:
            currentScanRes=currentScanRes/2  # make the scan finer in x
            startAbove    =(not startAbove)     # if we were above yvalue, switch to below and vice versa
            scanDirection =-1*scanDirection  # switch the x scan direction
            if abs(splineVal-yvalue) <= tolerance :
                upperLimit = scanX
                break

        lastPoint =splineVal
        scanX += scanDirection*currentScanRes

    scanX = startValue
    lastPoint=float('inf')
    currentScanRes=scanRes
    scanDirection=-1
    startAbove=spline.Eval(startValue) >= yvalue

    # Get the lower limit
    while (scanX > minValue) :
        splineVal=spline.Eval(scanX)

        # check if we've crossed the intersection yvalue
        switchCondition=(startAbove != (splineVal >= yvalue))

        if switchCondition:
            currentScanRes=currentScanRes/2  # make the scan finer in x
            startAbove    =(not startAbove)     # if we were above yvalue, switch to below and vice versa
            scanDirection =-1*scanDirection  # switch the x scan direction
            if abs(splineVal-yvalue) <= tolerance :
                lowerLimit = scanX
                break

        lastPoint = splineVal
        scanX += scanDirection*currentScanRes

    return lowerLimit,upperLimit


# now run and print
limList=options.limList.split(',')
fIn = ROOT.TFile("%s/%s"%(options.indir,options.iname),"READ")

tSpline=None
obsSpline=None

if options.splineV == "3" :
    tSpline=ROOT.TSpline3(fIn.Get("Spline%s%s"%(options.dist,options.prepost)))
    if options.unblind:
        obsSpline=ROOT.TSpline3(fIn.Get("SplineObs%s"%options.dist))
elif options.splineV == "2" :
    tSpline=ROOT.TMVA.TSpline2("Spline%s"%options.prepost, fIn.Get("CLsGraph%s%s"%(options.dist,options.prepost)))
    if options.unblind:
        obsSpline=ROOT.TMVA.TSpline2("SplineObs%s"%options.prepost, fIn.Get("CLsGraphincmlb"))



print "Limits for ",options.dist,' ',options.prepost

for lim in limList :
    lowLim,upLim=getSplineIntersection(float(lim),tSpline)
    print '\t\t(EXPECTED) For a CLs of ',lim,':  (',lowLim,' : ',upLim,')'

    if options.unblind :
        lowLim,upLim=getSplineIntersection(float(lim),obsSpline)
        print '\t\t(OBSERVED) For a CLs of ',lim,':  (',lowLim,' : ',upLim,')'



