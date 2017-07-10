import ROOT
import os
import sys
import optparse
import pickle
from collections import OrderedDict

from UESummaryPlotCommon import *

"""
"""
def buildUEPlot(obsAxis,sliceAxis,regions,fIn,fSyst):

    nSliceBins=1
    if sliceAxis: nSliceBins += sliceAxis.GetNbins()

    ueplot = UESummaryPlotInfo(obsAxis,sliceAxis)

    for i in xrange(0,nSliceBins):
      
        for r in regions:
    
            reg=str(r)

            #read the nominal expectations
            nomKey='%s_%s_%d_%s_None_True'%(ueplot.obs,ueplot.sliceVar,i,reg)
            t=fIn.Get(nomKey)
            print t,nomKey
            bkg=None
            for pkey in t.GetListOfKeys():
                h=t.Get(pkey.GetName())
                if not h.InheritsFrom('TH1') : continue
                if 'Data' in h.GetTitle(): 
                    ueplot.addData(h,i,r)
                elif h.GetTitle() in 't#bar{t}':
                    ueplot.addSignal(h,i,r)
                else:
                    if bkg is None: bkg=h.Clone('bkg')
                    else : bkg.Add(h)
                
            #subtract the background
            ueplot.subtractBackground(bkg,i,r)

#
#            #project experimental systematics and signal variations
#            signalVars.append( {} )
#            expSysts.append( {} )
#            expSystsKey='%s_%s_inc_syst_True'%(obs,s)
#            expSystsH=fIn.Get('%s/%s_%s'%(expSystsKey,expSystsKey,'t#bar{t}'))
#            for ybin in xrange(2,expSystsH.GetNbinsY()):
#                varName=expSystsH.GetYaxis().GetBinLabel(ybin)
#                systKey=varName
#                if systKey[-2:] in ['up','dn']  : systKey=systKey[:-2]
#                h=expSystsH.ProjectionX('px',ybin,ybin)
#                normalizePerSlice(h,obsAxis,sliceAxis)
#                if systKey in ['mur','muf','q'] :                    
#                    h.Divide(data)
#                    systKey='ME scale'
#                    if not systKey in signalVars[-1]: signalVars[-1][systKey]=[]
#                    signalVars[-1][systKey].append( h.Clone(varName) )
#                elif systKey in ['p_{T}(t)']:
#                    h.Divide(data)
#                    systKey='toppt'
#                    if not systKey in signalVars[-1]: signalVars[-1][systKey]=[]
#                    signalVars[-1][systKey].append( h.Clone(varName) )                    
#                elif systKey in ['tkeff','tkeffbcdef','tkeffgh','tkeffeta']:
#                    if systKey in ['tkeffbcdef','tkeffgh'] : continue
#                    h.Add(signal,-1)
#                    systKey='Trk. eff.'
#                    if not systKey in expSysts[-1]: expSysts[-1][systKey]=[]
#                    expSysts[-1][systKey].append( h.Clone(varName) )
#                else:
#                    h.Add(signal,-1)
#                    if not systKey in expSysts[-1]: expSysts[-1][systKey]=[]
#                    expSysts[-1][systKey].append(  h.Clone(varName) )
#
#            #variations to compare to
#            for varTitle in varTypes:
#                signalVars[-1][varTitle]=[]
#                for varName in varTypes[varTitle]:
#                    hvar=fSyst.Get(nomKey).Get(signal.GetName().replace('_nominal',' '+varName))
#                    normalizePerSlice(hvar,obsAxis,sliceAxis)
#                    hvar.Divide(data)
#                    signalVars[-1][varTitle].append(hvar)
#
    return ueplot

"""
"""
def readPlotsFrom(args,opt):

    outdir=args[0].replace('.root','')

    #analysis axes
    analysisaxis=None
    with open(opt.analysisAxis,'r') as cachefile:
        analysisaxis = pickle.load(cachefile)

    #open input files
    fIn=ROOT.TFile.Open(args[0])
    fSyst=ROOT.TFile.Open(args[1]) if len(args)>1 else None


    #build and show the different plot combinations (observables vs slices)
    for obs in OBSERVABLES:

        obsAxis=analysisaxis[(obs,True)]
        
        for s in SLICES:

            sliceAxis=None if s is None else analysisaxis[(s,False)]
            regions=['inc']
            if s in EVAXES: regions += [0,1,2]

            uePlot=buildUEPlot(obsAxis,sliceAxis,regions,fIn,fSyst)
            uePlot.show(outdir)

"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    #ROOT.gROOT.SetBatch(True)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--cfg',
                      dest='analysisAxis',
                      help='cfg with axis definitions [%default]',
                      type='string',
                      default='%s/src/TopLJets2015/TopAnalysis/UEanalysis/analysisaxiscfg.pck'%os.environ['CMSSW_BASE'])
    (opt, args) = parser.parse_args()

    readPlotsFrom(args,opt)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
