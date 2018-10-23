#!/usr/bin/env python

import ROOT
import optparse
import os,sys
from runDataFit import shushRooFit,showFitResult,observables,EVENTCATEGORIES

VAR2PLOT=[
    ('ly',  'Lepton rapidity',                  'y',                [-2.1, -1.0, -0.5, 0, 0.5, 1.0, 2.1]),
    ('lpt', 'Lepton transverse momentum [GeV]', 'p_{T} [GeV^{-1}]', [30,  40,  50,  60,  85, 200]),
    ]

def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    #ROOT.gROOT.SetBatch(True)
    shushRooFit()

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--output',    dest='output',    default='splot/',          type='string',   help='output directory [%default]')
    parser.add_option('-i', '--input',     dest='input',     default='fit_finalworkspace_wmodel0_0.root', type='string',   help='workspace [%default]')
    parser.add_option('-d', '--data',      dest='data',      default='pseudodata',      type='string',   help='data [%default]')
    parser.add_option('-v', '--verbose',   dest='verbose',   default=0,                            type=int,        help='Verbose mode [%default]')
    (opt, args) = parser.parse_args()

    #read workspace
    fIn=ROOT.TFile(opt.input)
    w=fIn.Get('w')

    #load the values obtained from the "standard" combined fit
    w.loadSnapshot('fitresult_combined')

    allVars=w.allVars()
    varIter = allVars.createIterator()
    var=varIter.Next()
    toFloat=['Nbkg_','xsec']    
    fixVarList=[]
    yieldsList=ROOT.RooArgList()
    while var :
        varName=var.GetName()    
        if not varName in ['mjj','mthad','mtlep','lpy','ly']: 
            fixThis=True
            for v in toFloat:
                if v in varName:
                    fixThis=False
                    break
            var.setConstant(fixThis)        
            if not fixThis:  
                yieldsList.add(var)
            else:
                fixVarList.append(varName)
        var=varIter.Next()
    print 'Floating this variables for the sPlot'
    yieldsList.Print()
    #print 'These variables have bin fixed',fixVarList

    #repeat the fit with a reduced set of variables
    inixsec=(w.var('xsec').getVal(),w.var('xsec').getError())
    pdf=w.pdf('model_combined_mjj')
    data=w.data(opt.data)

    #for cat in EVENTCATEGORIES:
    #    frame=w.var('mjj').frame()
    #    redData=data.reduce(ROOT.RooFit.Cut("sample==sample::%s"%cat))
    #    redData.plotOn(frame)
    #    pdf.plotOn(frame,ROOT.RooFit.ProjWData(redData))
    #    frame.Draw()
    #    raw_input(cat)


    pdf.fitTo(data,ROOT.RooFit.Extended())
   
    #nll=pdf.createNLL(data,ROOT.RooFit.Extended(True),ROOT.RooFit.NumCPU(8))
    #minuit=ROOT.RooMinuit(nll)
    #minuit.setStrategy(2)
    #minuit.migrad() #minimize with respect to all parameters
    #poi=ROOT.RooArgSet()    
    #poi.add(w.var('xsec'))
    #minuit.minos(poi)
    #r=minuit.save()

    #for cat in EVENTCATEGORIES:
    #    frame=w.var('mjj').frame()
    #    redData=data.reduce(ROOT.RooFit.Cut("sample==sample::%s"%cat))
    #    redData.plotOn(frame)
    #    pdf.plotOn(frame,ROOT.RooFit.ProjWData(redData))
    #    frame.Draw()
    #    raw_input(cat+' postfit')

    finalxsec=(w.var('xsec').getVal(),w.var('xsec').getError())
    print 'Initial fit gave the following cross section'
    print 'xsec=%3.3f +/- %3.3f'%inixsec
    print 'Reduced fit yields the following cross section'
    print 'xsec=%3.3f +/- %3.3f'%finalxsec
    
    #now do the splots (per category)
    iterator = yieldsList.createIterator()
    obj = iterator.Next()
    sigdata=None
    bkgperCat=[]
    while obj:
        objName=obj.GetName()
        if not 'xsec' in objName:
            catName=objName.replace('Nbkg_','')            
            print catName
            redData=data.reduce(ROOT.RooFit.Cut('sample==sample::%s'%catName))
            catYieldsList=ROOT.RooArgList(w.var('xsec'),w.var(objName))
            sData = ROOT.RooStats.SPlot("sdata","SPlotted data",redData,pdf,catYieldsList)
            print catName,sData.GetYieldFromSWeight('xsec_sw'),sData.GetYieldFromSWeight(objName+'_sw'),redData.sumEntries(),w.function(objName.replace('Nbkg','Nsig')).getVal()
            if sigdata:
                sigdata.append(ROOT.RooDataSet('signal'+catName,'signal',redData, redData.get(), '', 'xsec_sw'))
            else:
                sigdata=ROOT.RooDataSet('signal','signal',redData, redData.get(), '', 'xsec_sw')
                             
            bkgperCat.append( ROOT.RooDataSet('bkg'+catName, 'bkg '+catName,   redData, redData.get(), '', objName+'_sw') )
        obj = iterator.Next()
    sigdata.Print()
    
    print 'acc_mu',w.function('acc_mu').getVal()
    for v in ['eff_mu','lumi','xsec']:
        print v,w.var(v).getVal()
    
    #compare background vs signal like hypothesis
    os.system('mkdir -p %s'%opt.output)
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.02)
    c.SetBottomMargin(0.1)
    for var,varTitle,yvar,binDef in VAR2PLOT:        
        c.Clear()

        bins=ROOT.RooBinning(binDef[0],binDef[-1])
        for ibin in xrange(1,len(binDef)):
            bins.addBoundary(binDef[ibin])
        bins.Print('v')
        w.var(var).setBinning(bins)
        #w.var(var).setBins(5)
        frame=w.var(var).frame()
        for i in xrange(0,len(bkgperCat)):
            bkgperCat[i].plotOn(frame,
                                ROOT.RooFit.LineColor(ROOT.kGray),
                                ROOT.RooFit.LineStyle(1+i),
                                ROOT.RooFit.MarkerStyle(1+i),
                                ROOT.RooFit.MarkerColor(ROOT.kGray),
                                ROOT.RooFit.Name('bkg_%d'%i)
                                )
        sigdata.plotOn(frame,ROOT.RooFit.Name('sig'))
        frame.Draw()
        frame.GetYaxis().SetTitleOffset(1.3)
        frame.GetYaxis().SetTitle('dN/d%s'%yvar)
        frame.GetXaxis().SetTitle(varTitle)
        label = ROOT.TLatex()
        label.SetNDC()
        label.SetTextFont(42)
        label.SetTextSize(0.04)
        if opt.data=='pseudodata':
            label.DrawLatex(0.15,0.9,'#bf{CMS} #it{simulation preliminary}')
            label.DrawLatex(0.65,0.96,'#scale[0.8]{pPb (1pb^{-1}, #sqrt{s}=8.16 TeV)}')
        else:
            label.DrawLatex(0.15,0.9,'#bf{CMS} #it{preliminary}')
            label.DrawLatex(0.8,0.96,'(#scale[0.8]{#sqrt{s}=8.16 TeV})')
        #leg=ROOT.TLegend(0.15,0.88,0.4,0.75)
        #leg.SetFillColor(0)
        #leg.SetFillStyle(0)
        #leg.SetBorderSize(0)
        #leg.SetTextFont(42)
        #leg.SetTextSize(0.04)
        #leg.AddEntry('sig','Signal','ep')
        #for i in xrange(0,len(bkgperCat)):
        #    leg.AddEntry('bkg_%d'%i,bkperCat[i].GetTitle(),'l')
        #leg.Draw()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('%s/%s_splot.%s'%(opt.output,var,ext))



if __name__ == "__main__":
    main()
