#!/usr/bin/env python

import sys
import ROOT

#simple roofit version
def runFit(pdf,data,poi,constr=None):

    #create the log likelihood
    nll=pdf.createNLL(data,
                      ROOT.RooFit.Extended(True),
                      ROOT.RooFit.NumCPU(2))

    #add constraints
    if constr:
        parcels=ROOT.RooArgList(nll)
        iter = constr.createIterator()
        var = iter.Next()
        while var :
            parcels.add(var)
            var = iter.Next()
        nll=ROOT.RooAddition('nllc','nllc',parcels)

    minuit=ROOT.RooMinuit(nll)
    minuit.migrad() #minimize with respect to all parameters
    minuit.minos(poi)
    r=minuit.save()


fIn=ROOT.TFile.Open(sys.argv[1])

w=fIn.Get('w')
pdf=w.pdf('model_combined_mjj')
data=w.data('data')

poi=ROOT.RooArgSet()
poi.add(w.var('xsec'))

constr=ROOT.RooArgSet()
constr.add( w.function('ebconstraint') )

runFit(pdf,data,poi,constr)
raw_input()
