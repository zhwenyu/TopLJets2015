#!/usr/bin/env python

import ROOT
import optparse
import os,sys
from array import array
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

    #create a simplified PDF summing up all events with >=1 b
    catsOfInterest=['e1l4j1b1q','e1l4j2b','mu1l4j1b1q','mu1l4j2b']    
    data=w.data(opt.data).reduce(' || '.join(['sample==sample::%s'%c for c in catsOfInterest]))
    print 'Inclusive yields:',w.data(opt.data).sumEntries(),'->',data.sumEntries(),
    print 'after selecting only',catsOfInterest
    totalSig,totalBkg=0,0
    fSig,fBkg=[],[]
    sumSigExpr='SUM::Sgt1b('
    sumBkgExpr='SUM::Bgt1b('
    for i in xrange(0,len(catsOfInterest)):

        cat=catsOfInterest[i]

        nsig=w.function('Nsig_%s'%cat).getVal()
        totalSig+=nsig
        fSig.append(nsig)

        nbkg=w.var('Nbkg_%s'%cat).getVal()
        totalBkg+=nbkg
        fQCD=w.function('Nqcd_%s'%cat).getVal()/nbkg
        fBkg.append(nbkg*fQCD)
        fW=w.function('Nw_%s'%cat).getVal()/nbkg        
        fBkg.append(nbkg*fW)

        if i<len(catsOfInterest)-1:
            sumSigExpr += '{}*S_mjj_%s, '%cat
            sumBkgExpr += '{}*QCD_mjj_%s,{}*W_mjj_%s, '%('e' if cat[0]=='e' else 'mu',cat)
        else:
            sumSigExpr += 'S_mjj_%s'%cat
            sumBkgExpr += '{}*QCD_mjj_%s,W_mjj_%s'%('e' if cat[0]=='e' else 'mu',cat)

    fSig=tuple([x/totalSig for x in fSig][0:-1])
    fBkg=tuple([x/totalBkg for x in fBkg][0:-1])
    sumSigExpr += ' )'
    sumBkgExpr += ' )'
    sumSigExpr=sumSigExpr.format(*fSig)
    sumBkgExpr=sumBkgExpr.format(*fBkg)
    w.factory('nsiggt1b[%3.1f,%3.1f,%3.1f]'%(totalSig,totalSig*0.5,totalSig*2))
    w.factory(sumSigExpr)
    w.factory('nbkggt1b[%3.1f,%3.1f,%3.1f]'%(totalBkg,totalBkg*0.5,totalBkg*2))
    w.factory(sumBkgExpr)
    w.factory('SUM::modelgt1b( nsiggt1b*Sgt1b,nbkggt1b*Bgt1b )')    
    pdf=w.pdf('modelgt1b')

    print 'Simplified signal PDF:',sumSigExpr
    print 'Simplified background PDF:',sumBkgExpr    
    #frame=w.var('mjj').frame()
    #data.plotOn(frame)
    #pdf.plotOn(frame,ROOT.RooFit.ProjWData(data))
    #pdf.plotOn(frame,ROOT.RooFit.Components('Bgt1b'),ROOT.RooFit.LineColor(2)) ;
    #frame.Draw()

    allVars=w.allVars()
    varIter = allVars.createIterator()
    var=varIter.Next()
    toFloat=['nsiggt1b','nbkggt1b']    
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
    pdf.fitTo(data,ROOT.RooFit.Extended())
    
    iniyields=(totalSig,totalBkg)
    finalyields=(w.var('nsiggt1b').getVal(),w.var('nbkggt1b').getError())
    print 'Prefit S/B=(%3.1f/%3.1f)'%iniyields    
    print 'Postfit S/B=(%3.1f/%3.1f)'%finalyields
    

    #splot
    sData = ROOT.RooStats.SPlot("sdata","SPlotted data",data,pdf,yieldsList)
    sdset = sData.GetSDataSet()
    #sigdata=ROOT.RooDataSet('signal','signal',data, data.get(), '', 'nsiggt1b_sw')
    #bkgdata=ROOT.RooDataSet('bkg', 'bkg',data, data.get(), '', 'nbkggt1b_sw')
    
    #compare background vs signal like hypothesis
    os.system('mkdir -p %s'%opt.output)
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.1)
    dh=[]
    for var,varTitle,yvar,binDef in VAR2PLOT:        
        c.Clear()
        binArray = array('d',binDef)

        
        varList=ROOT.RooArgList(w.var(var))
        frame=w.var(var).frame()
        
        #for proc,marker in [('bkg',24),('sig',20)]:
        for proc,marker in [('sig',20)]:
            h=ROOT.TH1F(proc+var,proc+var,len(binArray)-1,binArray)
            for i in xrange(0,sdset.numEntries()):
                evargs=sdset.get(i)
                ival=evargs.find(var).getVal()
                sw=sData.GetSWeight(i,'n%sgt1b'%proc)
                xbin=h.GetXaxis().FindBin(ival)
                wgt=sw/h.GetXaxis().GetBinWidth(xbin)
                h.Fill(ival,wgt)
            dh.append( ROOT.RooDataHist(proc+var+'dh',proc+var+'dh',varList,h) )
            dh[-1].plotOn(frame,ROOT.RooFit.MarkerStyle(marker),ROOT.RooFit.Name(proc))
            h.Delete()

        frame.Draw()
        frame.GetYaxis().SetTitleSize(0.04)
        frame.GetXaxis().SetTitleSize(0.04)
        frame.GetYaxis().SetTitleOffset(1.5)
        frame.GetYaxis().SetTitle('dN/d%s'%yvar)
        frame.GetXaxis().SetTitle(varTitle)
        frame.GetXaxis().SetRangeUser(binDef[0],binDef[-1])
        frame.GetYaxis().SetRangeUser(0,frame.GetMaximum()*1.3)
        label = ROOT.TLatex()
        label.SetNDC()
        label.SetTextFont(42)
        label.SetTextSize(0.04)
        if opt.data=='pseudodata':
            label.DrawLatex(0.17,0.9,'#bf{CMS} #it{simulation preliminary}')
            label.SetTextAlign(32)
            label.DrawLatex(0.95,0.97,'#scale[0.8]{pPb (2pb^{-1}, #sqrt{s}=8.16 TeV)}')
        else:
            label.DrawLatex(0.15,0.9,'#bf{CMS} #it{preliminary}')
            label.SetTextAlign(32)
            label.DrawLatex(0.95,0.97,'(#scale[0.8]{#sqrt{s}=8.16 TeV})')
        """
        leg=ROOT.TLegend(0.17,0.88,0.4,0.75)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        leg.AddEntry('sig','t#bar{t}','ep')
        leg.AddEntry('bkg','total background','ep')
        leg.Draw()
        """
        c.Modified()
        c.Update()
        for ext in ['png','pdf','root']:
            c.SaveAs('%s/%s_splot.%s'%(opt.output,var,ext))

    fOut=ROOT.TFile.Open('splots.root','RECREATE')
    dh[0].createHistogram("splot_ly",w.var('ly')).Write()
    dh[1].createHistogram("splot_lpt",w.var('lpt')).Write()
    fOut.Close()

if __name__ == "__main__":
    main()

