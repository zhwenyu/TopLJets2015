#!/usr/bin/env python2.7

#from prepareWorkspace import EVENTCATEGORIES
EVENTCATEGORIES=['1l4j2b','1l4j2q','1l4j1b1q']

import ROOT
import optparse
import os,sys
from roofitTools import *
from parameterizeMCShapes import ALLPDFS

"""
disable RooFit verbosity
"""
def shushRooFit():
    ROOT.RooMsgService.instance().setSilentMode(True);
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Minimization)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Minimization)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.DataHandling)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Fitting)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Plotting)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Plotting)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.InputArguments)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.InputArguments)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Eval)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Eval)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Integration)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)

"""
"""
def addToFitResults(title,fitResults,varSet):
    fitResults[title]={}
    iter = varSet.createIterator()
    iparam = iter.Next()
    while iparam :
        if not iparam.getAttribute('Constant'):
            fitResults[title][iparam.GetName()]=(iparam.getVal(),iparam.getError())
        iparam = iter.Next()

"""
"""
def printFitResults(fitResults,opt):
    allVars=set()
    for key in fitResults:
        for pname in fitResults[key]:
            allVars.add(pname)

    with open('%s/fitResults_%d.tex'%(opt.output,opt.onlyResonant),'w') as fOut:
        fOut.write('\\hline\n')
        fOut.write('\multirow{3}{*}{Variable} & \\multicolumn{3}{c}{Fit type} \\\\\n')
        for key in fitResults:
            fOut.write('& %s'%key)
        fOut.write('\\\\\n')
        fOut.write('\\hline\n')
        for pname in allVars:
            fOut.write('%25s'%pname)
            for key in fitResults:
                if not pname in fitResults[key]: 
                    fOut.write('& %25s'%'')
                else : 
                    fOut.write('& %25s'%'$%3.2f \\pm %3.2f$'%fitResults[key][pname])
            fOut.write('\\\\\n')
        fOut.write('\\hline\n')

"""
"""
def defineCombinedPDFs(w):
    for t in ['2D','3D']:
        simPDF=ROOT.RooSimultaneous('model_%s'%t,'model_%s'%t,w.cat('sample'))
        for cat in EVENTCATEGORIES:

            #multiply PDFs
            sigProd='S_mjj_{0},S_mthad_{0}'.format(cat)
            qcdProd='QCD_mjj,QCD_mthad'
            wProd  ='W_mjj_{0},W_mthad_{0}'.format(cat)
            if t=='3D':
                sigProd += ',S_mtlep_{0}'.format(cat)
                qcdProd += ',QCD_mtlep'
                wProd   += ',W_mtlep_{0}'.format(cat)
            w.factory('PROD::S_{1}_{0}({2})'.format(cat,t,sigProd))
            w.factory('PROD::QCD_{0}({1})'.format(t,qcdProd))
            w.factory('PROD::W_{1}_{0}({2})'.format(cat,t,wProd))

            #sum up so that it can be extended
            w.factory('SUM::model_{1}_{0}(Nsig_{0}*S_{1}_{0},Nqcd_{0}*QCD_{1},Nw_{0}*W_{1}_{0})'.format(cat,t))
            simPDF.addPdf(w.pdf('model_{1}_{0}'.format(cat,t)),cat)

        #import simultaneous pdf to workspace
        getattr(w,'import')(simPDF,ROOT.RooFit.RecycleConflictNodes())

"""
"""
def definePDF(w,varName):

    # PDF definition
    simPDF=ROOT.RooSimultaneous('model_%s'%varName,'model_%s'%varName,w.cat('sample'))

    #YIELDS
    #signal yields are parametrized as Nsig=fsig2.Nsig+ fsig1.Nsig + (1-fsig1-fsig2).Nsig
    w.factory("Nsig[500,0,5000]")
    #w.factory("RooFormulaVar::Nsig_1l4j2b('@0*@1',{fsig2[0,1],Nsig})")
    #w.factory("RooFormulaVar::Nsig_1l4j1b1q('@0*@1',{fsig1[0,1],Nsig})")
    #w.factory("RooFormulaVar::Nsig_1l4j2q('(1-@0-@1)*@2',{fsig2,fsig1,Nsig})")
    w.factory("RooFormulaVar::Nsig_1l4j2b('@0*pow(@1,2)',{Nsig,eb[0.65,0.9]})")
    w.factory("RooFormulaVar::Nsig_1l4j1b1q('@0*2*@1*(1-@1)',{Nsig,eb})")
    w.factory("RooFormulaVar::Nsig_1l4j2q('@0*pow(1-@1,2)',{Nsig,eb})")

    #backgrounds yields are indepedent in each category (to be profiled)
    #the number of QCD (W) like is parameterized as fqcd.Nbkg ((1-fqcd).Nw)
    for cat in EVENTCATEGORIES:

        w.factory('Nbkg_{0}[1000,0,100000]'.format(cat))
        if '2b'   in cat : w.factory('fqcd_{0}[0.0,0.9]'.format(cat))
        elif '1b' in cat : w.factory('fqcd_{0}[0.0,0.9]'.format(cat))
        else             : w.factory('fqcd_{0}[0.0,0.9]'.format(cat))

        w.factory("RooFormulaVar::Nqcd_{0}('@0*@1',{{Nbkg_{0},fqcd_{0}}})".format(cat))

        w.factory("RooFormulaVar::Nw_{0}('@0*(1-@1)',{{Nbkg_{0},fqcd_{0}}})".format(cat))
        minMPV=10 if varName=='mjj' else 150
        maxMPV=60 if varName=='mjj' else 200
        w.factory('RooLandau::W_{1}_{0}({1},mpv_w{1}_{0}[{2},{3}],width_w{1}_{0}[10,100])'.format(cat,varName,minMPV,maxMPV))

        #w.factory('SUM:model_{0}(Nqcd_{0}*QCD_{1},Nw_{0}*W_{0})'.format(cat,varName))
        w.factory('SUM::model_{1}_{0}(Nsig_{0}*S_{1}_{0},Nqcd_{0}*QCD_{1},Nw_{0}*W_{1}_{0})'.format(cat,varName))
        simPDF.addPdf(w.pdf('model_{1}_{0}'.format(cat,varName)),cat)

    #import simultaneous pdf to workspace
    getattr(w,'import')(simPDF,ROOT.RooFit.RecycleConflictNodes())
    return w.pdf('model_%s'%varName)


"""
model QCD using RooKeysPDF and the failing id/iso/btag sideband data
"""
def addQCDModel(w,opt,varName):

    #get sideband data
    redData = w.data('data').reduce(ROOT.RooFit.Cut("sample==sample::1f4j2q"))
    fullsideBand=redData
    if redData.numEntries()>5000:
        print 'reducing numentries to 5000'
        redData=redData.reduce(ROOT.RooFit.EventRange(0,5000))

    #run RooKeysPdf and add it to the workspace
    keyspdf=ROOT.RooKeysPdf("QCD_%s"%varName,"QCD_%s"%varName,
                            w.var(varName),
                            redData,
                            ROOT.RooKeysPdf.NoMirror,
                            1.4)
    getattr(w,'import')(keyspdf,ROOT.RooFit.RecycleConflictNodes())

    #save fit result plot
    showFitResult(fitVar=varName,
                  data=fullsideBand,
                  pdf=w.pdf('QCD_%s'%varName),
                  categs=[''],
                  w=w,
                  showComponents=[],
                  rangeX=(0,400),
                  tagTitle='QCD',
                  outDir=opt.output)

"""
readout signal models
"""
def addSignalModel(w,opt,varName):

    #get signal workspace from file
    fIn=ROOT.TFile.Open(opt.signal)
    ws=fIn.Get('w')
    fIn.Close()

    #define resolution
    #w.factory('bias_s[0]')
    #w.var('bias_s').SetTitle('#mu_{bias}')
    #w.factory('resol_s[0.1]')
    #w.var('resol_s').SetTitle('#sigma_{resol}')
    #w.factory("Gaussian::S_resolution(%s,bias_s,resol_s)"%opt.varName)

    data=w.data('data')
    for cat in EVENTCATEGORIES:

        if opt.onlyResonant<2:
            for c in ['cor','wro']:
                pdfDef=[x.format(cat) for x in  ALLPDFS[('S_%s'%c,varName)] ]
                for p in pdfDef: pdf=w.factory(p)
                iter = pdf.getParameters(data).createIterator()
                iparam = iter.Next()
                while iparam :
                    mcval=ws.var(iparam.GetName()).getVal()
                    iparam.setVal(mcval)
                    iparam.setConstant(True)
                    iparam = iter.Next()

            resFracName='f_s{1}_{0}'.format(cat,varName)
            w.factory('SUM::S_{1}_{0}( {2}[0,1]*S_cor{1}_{0}, S_wro{1}_{0} )'.format(cat,varName,resFracName))
            w.var(resFracName).setVal(1.0 if opt.onlyResonant==1 else ws.var(resFracName).getVal() )
            w.var(resFracName).setConstant(True)
        else:
            muLimits=[60,90] if varName=='mjj' else [150,180]
            sigmaLimits=[10,30]
            w.factory('mu_scor%s[60,90]'%varName)
            w.factory('sigmaL_scor%s[10,25]'%varName)
            w.factory('sigmaR_scor%s[10,25]'%varName)
            w.factory('RooBifurGauss::S_{1}_{0}({1},mu_scor{1},sigmaL_scor{1},sigmaR_scor{1})'.format(cat,varName))

        #convolve with a common resolution
        #w.factory("FCONV::S_{0}({1}, S_resolution,S_base_{0})".format(cat,opt.varName))

        #frame=w.var(opt.varName).frame()
        #pdf.plotOn(frame)
        #w.pdf('S_{0}'.format(cat)).plotOn(frame)
        #frame.Draw()
        #raw_input()



"""
minimizes the likelihood and profiles the poi
"""
def runFit(pdf,data,poi,obs,w,outDir):

    mc=ROOT.RooStats.ModelConfig('mc_%s'%pdf.GetName(),w)
    mc.SetPdf(pdf)
    mc.SetParametersOfInterest(poi)
    mc.SetObservables(obs)
    getattr(w,'import')(mc)

    pl=ROOT.RooStats.ProfileLikelihoodCalculator(data,mc)
    pl.SetConfidenceLevel(0.683);
    interval=pl.GetInterval()

    ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")
    firstPOI = mc.GetParametersOfInterest().first()
    lowerLimit = interval.LowerLimit(firstPOI)
    upperLimit = interval.UpperLimit(firstPOI)
    cenVal     = interval.GetBestFitParameters().find(firstPOI.GetName()).getVal()
    unc=(upperLimit-lowerLimit)*0.5

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.1)
    plot = ROOT.RooStats.LikelihoodIntervalPlot(interval)
    plot.SetRange(lowerLimit-unc*2,upperLimit+unc*2)
    plot.SetNPoints(40)
    plot.Draw("tf1")
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.04)
    label.DrawLatex(0.6,0.9,'#bf{CMS} #it{preliminary}')
    label.DrawLatex(0.6,0.85,'%s=%3.1f^{+%3.1f}_{-%3.1f}'%(firstPOI.GetTitle(),cenVal,upperLimit-cenVal,cenVal-lowerLimit))
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s/ll_%s.%s'%(outDir,mc.GetName(),ext))

#simple roofit version
def runSimpleFit(pdf,data,poi):
    nll=pdf.createNLL(data,ROOT.RooFit.Extended(True),ROOT.RooFit.NumCPU(2))
    ROOT.RooMinimizer(nll).migrad() #minimize with respect to all parameters
    pll=nll.createProfile(poi)
    return (nll,pll)

"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(False) #True)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--output',    dest='output',    default='plots/Data8TeV_pp',          type='string',   help='output directory [%default]')
    parser.add_option('-i', '--input',     dest='input',     default='workspace_Data8TeV_pp.root', type='string',   help='workspace [%default]')
    parser.add_option('-s', '--signal',    dest='signal',    default='pdf_workspace_MC8.16TeV_TTbar_pPb.root', type='string',   help='signal workspace [%default]')
    parser.add_option(      '--onlyResonant',   dest='onlyResonant',   default=0,                  type=int,        help='0-full signal;1-res from MC;2-res from CB [%default]')
    parser.add_option('-v', '--verbose',   dest='verbose',   default=0,                            type=int,        help='Verbose mode [%default]')
    (opt, args) = parser.parse_args()

    if opt.verbose<9 : shushRooFit()

    #read workspace
    fIn=ROOT.TFile(opt.input)
    w=fIn.Get('w')
    fIn.Close()
    observables=[('mjj','M(jj)'),('mthad','M(t_{had})'),('mtlep','M(t_{lep})')]
    for vname,vtit in observables:
        w.var(vname).SetTitle(vtit)

        #model QCD
        addQCDModel(w,opt,vname)

        #add signal model
        addSignalModel(w,opt,vname)

        #define the PDFs
        definePDF(w,vname)

    fitResults={}

    #parameter of interest
    w.var('Nsig').SetTitle('N(S)')
    poi=ROOT.RooArgSet()
    poi.add(w.var('Nsig'))

    #data to fit
    data=w.data('data')

    #run the fit to the mjj variable
    pdf=w.pdf('model_mjj')
    obsList=ROOT.RooArgSet()
    obsList.add(w.var('mjj'))
    #runFit(pdf,data,poi,obsList,w,opt.output)
    runSimpleFit(pdf,data,poi)
    addToFitResults('1D',fitResults,pdf.getParameters(data))

    compsToShow=['S_cor*','S_cor*,S_wro*','S_cor*,S_wro*,W_*']
    if opt.onlyResonant==1 : compsToShow=['S_cor*','S_cor*,W_*']
    if opt.onlyResonant==2 : compsToShow=['S_*','S_*,W_*']
    paramList=[('Nsig','N(S)')] #,('eb','#varepsilon_{b}')]
    if opt.onlyResonant==2 : paramList += [ ('mu_scormjj','#mu'),
                                            ('sigmaL_scormjj','#sigma_{L}'),
                                            ('sigmaR_scormjj','#sigma_{R}')]
    showFitResult(fitVar='mjj',
                  data=data,
                  pdf=pdf,
                  categs=EVENTCATEGORIES,
                  w=w,
                  showComponents=compsToShow,
                  rangeX=(0,300),
                  outDir=opt.output,
                  paramList=paramList,
                  pfix='_%dfit'%opt.onlyResonant)

    #run the 2D/3D fits
    if opt.onlyResonant==0:

        defineCombinedPDFs(w)
        for t in ['2D','3D']:
            pdf=w.pdf('model_{0}'.format(t))
            obsList=ROOT.RooArgSet()
            obsList.add(w.var('mjj'))
            obsList.add(w.var('mthad'))
            if t=='3D': obsList.add(w.var('mtlep'))
            #runFit(pdf,data,poi,obsList,w,opt.output)
            runSimpleFit(pdf,data,poi)
            addToFitResults(t,fitResults,pdf.getParameters(data))
            for obs in observables:
                if obs=='mtlep' and t!='3D' : continue
                showFitResult(fitVar=obs[0],
                              data=data,
                              pdf=pdf,
                              categs=EVENTCATEGORIES,
                              w=w,
                              showComponents=compsToShow,
                              rangeX=(0,300 if obs[0]=='mjj' else 400),
                              outDir=opt.output,
                              paramList=paramList,
                              pfix='_%sfit'%t)

    printFitResults(fitResults,opt)


if __name__ == "__main__":
    main()
