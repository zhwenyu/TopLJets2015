#!/usr/bin/env python

import ROOT
import optparse
from roofitTools import *
from prepareWorkspace import EVENTCATEGORIES as SELEVENTCATEGORIES
EVENTCATEGORIES=[x for x in SELEVENTCATEGORIES if not '1f' in x]


"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True) #False)
    ROOT.gErrorIgnoreLevel = 5000

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',     dest='input',     default='finalfitworkskace_0.root', type='string',   help='workspace [%default]')
    parser.add_option('-d', '--data',      dest='data',      default='data', type='string',   help='data name [%default]')
    (opt, args) = parser.parse_args()

    fIn=ROOT.TFile.Open(opt.input)
    w=fIn.Get('w')
    w.var('mjj').SetTitle("m_{jj'} [GeV]")
    w.var('mthad').SetTitle('m_{top} [GeV]')
    w.var('mtlep').SetTitle('m_{l#nub} [GeV]')

    data=w.data(opt.data)

    #fix all parameters related to mjj
    w.loadSnapshot('fitresult_combined')
    iter = w.pdf('model_combined_mjj').getParameters(data).createIterator()
    iparam = iter.Next()
    while iparam :
        iparam.setConstant(True)
        iparam = iter.Next()

    #fix background parameters for mthad and mtlep
    #pdf3d=w.pdf('model_combined_3D')
    #pdf3d.fitTo(data)
    pdfmthad=w.pdf('model_combined_mthad')
    pdfmthad.fitTo(data)
    pdfmtlep=w.pdf('model_combined_mtlep')
    pdfmtlep.fitTo(data)

    #compsToShow=[('S_cor*','t#bar{t} cor. perm.'),('S_cor*,S_wro*','t#bar{t} wro. perm.')] #,'S_cor*,S_wro*,W_*']
    compsToShow=[('S_wro*,S_cor2*','t#bar{t} wrong'),('S_cor*,S_wro*','t#bar{t} correct')] #,'S_cor*,S_wro*,W_*']
    compsToShow=[('QCD_*,W_*','background'),('QCD_*,W_*,S_wro*,S_cor2*','t#bar{t} wrong')] #,'S_cor*,S_wro*,W_*']
    
    results={}
    for x in EVENTCATEGORIES:

        xtitle  = 'e' if x[0]=='e' else '#mu'
        xtitle += '^{#pm} + #geq4j'
        if '2q' in x   : xtitle += ' (=0b)'
        if '1b1q' in x : xtitle += ' (=1b)'
        if '2b' in x   : xtitle += ' (#geq2b)'
        #xtitle += ' events'

        for var,rangeX,model in [('mjj',(25,300),'model_combined_mjj'),                                 
                                 #('mthad',(100,400),'model_combined_mthad'),
                                 #('mtlep',(100,400),'model_combined_mtlep')
                                 ]:
            if var=='mjj' : w.var(var).SetTitle("m_{jj'} [GeV]")
            for noPulls in [True, False]:
                res=showFitResult(fitVar=var,
                                  data=data,
                                  pdf=w.pdf(model),
                                  categs=[x],
                                  tagTitle=xtitle,
                                  w=w,
                                  showComponents=compsToShow,
                                  rangeX=rangeX,
                                  outDir='./',
                                  paramList=None,
                                  pfix='_final',
                                  extsToSave=['png','pdf','root'],
                                  noPulls=noPulls)
            results[(res[0],res[1])]=res[2]
            #raw_input()

    import pickle
    pickle.dump( results, open( "chisquare_final.pck", "wb" ) )
        

if __name__ == "__main__":
    main()
