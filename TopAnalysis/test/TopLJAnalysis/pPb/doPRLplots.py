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
    (opt, args) = parser.parse_args()

    fIn=ROOT.TFile.Open(opt.input)
    w=fIn.Get('w')
    w.var('mjj').SetTitle('M(jj) [GeV]')
    w.var('mthad').SetTitle('M_{top} [GeV]')
    w.var('mtlep').SetTitle('M(t_{lep}) [GeV]')

    #fix all parameters related to mjj
    w.loadSnapshot('fitresult_combined')
    iter = w.pdf('model_combined_mjj').getParameters(w.data('data')).createIterator()
    iparam = iter.Next()
    while iparam :
        print iparam.GetName()
        iparam.setConstant(True)
        iparam = iter.Next()

    #fix background parameters for mthad and mtlep
    #pdf3d=w.pdf('model_combined_3D')
    #pdf3d.fitTo(w.data('data'))
    pdfmthad=w.pdf('model_combined_mthad')
    pdfmthad.fitTo(w.data('data'))
    pdfmtlep=w.pdf('model_combined_mtlep')
    pdfmtlep.fitTo(w.data('data'))

    #compsToShow=[('S_cor*','t#bar{t} cor. perm.'),('S_cor*,S_wro*','t#bar{t} wro. perm.')] #,'S_cor*,S_wro*,W_*']
    compsToShow=[('S_wro*','t#bar{t} wrong perm.'),('S_cor*,S_wro*','t#bar{t} correct perm.')] #,'S_cor*,S_wro*,W_*']
    
    for x in EVENTCATEGORIES:

        xtitle  = 'e' if x[0]=='e' else '#mu'
        xtitle += '^{#pm} + #geq4j'
        if '2q' in x   : xtitle += ' (=0b)'
        if '1b1q' in x : xtitle += ' (=1b)'
        if '2b' in x   : xtitle += ' (#geq2b)'
        xtitle += ' events'

        for var,rangeX,model in [('mjj',(25,300),'model_combined_mjj'),                                 
                                 ('mthad',(100,400),'model_combined_mthad'),
                                 ('mtlep',(100,400),'model_combined_mtlep')
                                 ]:
            showFitResult(fitVar=var,
                        data=w.data('data'),
                        pdf=w.pdf(model),
                        categs=[x],
                        tagTitle=xtitle,
                        w=w,
                        showComponents=compsToShow,
                        rangeX=rangeX,
                        outDir='./',
                        paramList=None,
                        pfix='_final',
                        extsToSave=['png','pdf','root'])


        

if __name__ == "__main__":
    main()
