import sys
import ROOT
from compareOptimResults import showShapes

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

bosonName=sys.argv[3]
if bosonName=='g' : bosonName='#gamma'

showShapes(resultsDir=sys.argv[1],
           name='shapes',
           title='pp%sX'%bosonName,
           mass=int(sys.argv[2]),
           boson=sys.argv[3],
           r95=None,
           sig=None,
           lumi=37.500,
           showPseudoData=False,
           showAllBkgs=True) #False)
