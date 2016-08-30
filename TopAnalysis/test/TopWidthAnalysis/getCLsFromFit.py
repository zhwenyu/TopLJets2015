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
    epilog="Collects confidence level information from CLs fits"
    )
parser.add_option("-i",     type="string", dest="indir"  , default="./",     help="directory to look for stats files in")
parser.add_option("--iname",type="string", dest="iname"  , default="statsPlots.root", help="filename containing fit info")
parser.add_option("--dist", type="string", dest="dist"   , default="incmlb",    help="the observable distribution to look at")
parser.add_option("--nomWid", type="string", dest="inwid"   , default="1.324",    help="the observable distribution to look at")
parser.add_option("--unblind", dest="unblind", default=False, action='store_true',  help="the observable distribution to look at")

(options, args) = parser.parse_args()


def getSplineIntersection(yvalue, spline, scanRes=0.01, startValue=float(options.inwid),minValue=0,maxValue=6) :

    scanX=startValue
    lastDistance=float('inf')
    upperLimit=None
    lowerLimit=None

    # Get the upper limit
    while (scanX < maxValue) :
        distToSpline=spline.DistancetoPrimitive(scanX,yvalue)

        if distToSpline >= lastDistance:
            upperLimit = scanX - scanRes/2
            break

        lastDistance = distToSpline
        scanX += scanRes

    scanX = startValue
    lastDistance=float('inf')

    # Get the lower limit
    while (scanX > minValue) :
        distToSpline=spline.DistancetoPrimitive(scanX,yvalue)

        if distToSpline >= lastDistance:
            lowerLimit = scanX + scanRes/2
            break

        lastDistance = distToSpline
        scanX += -1*scanRes

    return lowerLimit,upperLimit


# now run and print
limList=options.limList.split(',')
fIn = ROOT.TFile("%s/%s"%(options.indir,options.iname),"READ"
tSpline=fIn.Get("Splinepost%s"%options.dist)
obsSpline=None

if options.unblind:
    obsSpline=fIn.Get("SplineObs%s"%options.dist)

for lim in limList :
    lowLim,upLim=getSplineIntersection(float(lim),tSpline)
    print '\t\t(EXPECTED) For a confidence level of ',lim,':  (',lowLim,' : ',upLim,')'

    if options.unblind :
        lowLim,upLim=getSplineIntersection(float(lim),obsSpline)
        print '\t\t(OBSERVED) For a confidence level of ',lim,':  (',lowLim,' : ',upLim,')'



