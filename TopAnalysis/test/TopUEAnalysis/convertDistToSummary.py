import ROOT
import sys
import pickle
import os
from UEAnalysisHandler import VARTITLES
from UEPlot import *

analysisCfg=sys.argv[1]
with open(analysisCfg,'r') as cachefile:
    obsAxis = pickle.load(cachefile)[('gen','axis')]
    cuts    = pickle.load(cachefile)
    obs     = pickle.load(cachefile)

inROOT=sys.argv[2]
fIn=ROOT.TFile(inROOT)
uePlots={}
uePlots['PW+PY8']=UEPlot(obs,VARTITLES[obs],obsAxis)
h=fIn.Get('gen')
uePlots['PW+PY8'].addVariation('PW+PY8',None,h)
uePlots['PW+PY8'].finalize()

with open(os.path.join( os.path.dirname(analysisCfg),'unfold/unfold_summary.pck'),'w') as cachefile:
    pickle.dump( uePlots, cachefile, pickle.HIGHEST_PROTOCOL)
