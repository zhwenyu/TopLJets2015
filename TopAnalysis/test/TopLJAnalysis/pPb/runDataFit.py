#!/usr/bin/env python2.7

from prepareWorkspace import EVENTCATEGORIES as SELEVENTCATEGORIES
EVENTCATEGORIES=[x for x in SELEVENTCATEGORIES if not '1f' in x]

lumi=174.5
acceptance={'e':0.056,'mu':0.060}
efficiency={'e':0.95*0.85*0.95,'mu':0.954*0.993*0.985*0.981}
ebExp,ebUnc=0.595,0.10*0.595


import ROOT
import optparse
import os,sys
from roofitTools import *
from parameterizeMCShapes import ALLPDFS

observables=[('mjj','M(jj)'),('mthad','M(t_{had})'),('mtlep','M(t_{lep})')]

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

    with open('%s/fitResults_%d.tex'%(opt.output,opt.fitType),'w') as fOut:
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

        simPDF={}
        for ch in ['combined','e','mu']:
            simPDF[ch]=ROOT.RooSimultaneous('model_%s_%s'%(ch,t),'model_%s_%s'%(ch,t),w.cat('sample'))

        for cat in EVENTCATEGORIES:

            ch='e' if cat[0]=='e' else 'mu'

            #multiply PDFs
            sigProd='S_mjj_{0},S_mthad_{0}'.format(cat)
            qcdProd='QCD_mjj_{0},QCD_mthad_{0}'.format(ch)
            wProd  ='W_mjj_{0},W_mthad_{0}'.format(cat)
            if t=='3D':
                sigProd += ',S_mtlep_{0}'.format(cat)
                qcdProd += ',QCD_mtlep_{0}'.format(ch)
                wProd   += ',W_mtlep_{0}'.format(cat)
            w.factory('PROD::S_{1}_{0}({2})'.format(cat,t,sigProd))
            w.factory('PROD::QCD_{0}_{1}({2})'.format(t,ch,qcdProd))
            w.factory('PROD::W_{1}_{0}({2})'.format(cat,t,wProd))

            #sum up so that it can be extended
            w.factory('SUM::model_{1}_{0}(Nsig_{0}*S_{1}_{0},Nqcd_{0}*QCD_{1}_{2},Nw_{0}*W_{1}_{0})'.format(cat,t,ch))
            simPDF['combined'].addPdf(w.pdf('model_{1}_{0}'.format(cat,t)),cat)
            simPDF[ch].addPdf(w.pdf('model_{1}_{0}'.format(cat,t)),cat)

        #import simultaneous pdf to workspace
        for ch in ['combined','e','mu']:
            print 'Importing',simPDF[ch].GetName()
            getattr(w,'import')(simPDF[ch],ROOT.RooFit.RecycleConflictNodes())

"""
"""
def definePDF(w,varName):

    # PDF definition
    simPDF={'combined':ROOT.RooSimultaneous('model_combined_%s'%varName,'model_combined_%s'%varName,w.cat('sample'))}

    #YIELDS
    for ch in ['e','mu']:
        simPDF[ch]=ROOT.RooSimultaneous('model_%s_%s'%(ch,varName),'model_%s_%s'%(ch,varName),w.cat('sample'))
        w.factory("RooFormulaVar::Nsig_{0}('@0*@1*@2*@3',{{acc_{0},eff_{0},lumi,xsec}})".format(ch))
        w.factory("RooFormulaVar::Nsig_{0}1l4j2b('@0*pow(@1,2)',{{Nsig_{0},eb}})".format(ch))
        w.factory("RooFormulaVar::Nsig_{0}1l4j1b1q('@0*2*@1*(1-@1)',{{Nsig_{0},eb}})".format(ch))
        w.factory("RooFormulaVar::Nsig_{0}1l4j2q('@0*pow(1-@1,2)',{{Nsig_{0},eb}})".format(ch))

    #backgrounds yields are indepedent in each category (to be profiled)
    #the number of QCD (W) like is parameterized as fqcd.Nbkg ((1-fqcd).Nw)
    for cat in EVENTCATEGORIES:

        w.factory('Nbkg_{0}[1000,0,100000]'.format(cat))
        if '2b'   in cat : w.factory('fqcd_{0}[0.0,0.5]'.format(cat))
        elif '1b' in cat : w.factory('fqcd_{0}[0.0,0.9]'.format(cat))
        else             : w.factory('fqcd_{0}[0.0,1.0]'.format(cat))
        #if '2b'   in cat : w.factory('fqcd_{0}[0.0]'.format(cat))
        #elif '1b' in cat : w.factory('fqcd_{0}[0.0]'.format(cat))
        #else             : w.factory('fqcd_{0}[0.0]'.format(cat))

        w.factory("RooFormulaVar::Nqcd_{0}('@0*@1',{{Nbkg_{0},fqcd_{0}}})".format(cat))

        #W background shares common parameters for shape in electron and muon samples
        w.factory("RooFormulaVar::Nw_{0}('@0*(1-@1)',{{Nbkg_{0},fqcd_{0}}})".format(cat))
        ch='e' if cat[0]=='e' else 'mu'
        baseCat=cat[1:] if cat[0]=='e' else cat[2:]
        mpvName,widthName='mpv_w{1}_{0}'.format(baseCat,varName),'width_w{1}_{0}'.format(baseCat,varName)
        minMPV=20  if varName=='mjj' else 150
        maxMPV=80 if varName=='mjj' else 200
        w.factory('{0}[{1},{2}]'.format(mpvName,minMPV,maxMPV))
        w.factory('{0}[10,100]'.format(widthName))
        w.factory('RooLandau::W_{1}_{0}({1},{2},{3})'.format(cat,varName,mpvName,widthName))

        #w.factory('SUM:model_{0}(Nqcd_{0}*QCD_{1},Nw_{0}*W_{0})'.format(cat,varName))
        w.factory('SUM::model_{1}_{0}(Nsig_{0}*S_{1}_{0},Nqcd_{0}*QCD_{1}_{2},Nw_{0}*W_{1}_{0})'.format(cat,varName,ch))
        w.cat('sample').setLabel(cat)
        simPDF['combined'].addPdf(w.pdf('model_{1}_{0}'.format(cat,varName)),cat)
        simPDF[ch].addPdf(w.pdf('model_{1}_{0}'.format(cat,varName)),cat)


    #import simultaneous pdf to workspace
    for ch in simPDF:
        print 'Importing',simPDF[ch].GetName()
        getattr(w,'import')(simPDF[ch],ROOT.RooFit.RecycleConflictNodes())

"""
model QCD using RooKeysPDF and the failing id/iso/btag sideband data
"""
def addQCDModel(w,opt,varName):

    #get sideband data
    for ch in ['e','mu']:
        redData = w.data('data').reduce(ROOT.RooFit.Cut("sample==sample::%s1f4j2q"%ch))
        fullsideBand=redData
        if redData.numEntries()>5000:
            print 'reducing numentries to 5000'
            redData=redData.reduce(ROOT.RooFit.EventRange(0,5000))

        #run RooKeysPdf and add it to the workspace
        keyspdf=ROOT.RooKeysPdf("QCD_%s_%s"%(varName,ch),
                                "QCD_%s_%s"%(varName,ch),
                                w.var(varName),
                                redData,
                                ROOT.RooKeysPdf.NoMirror,
                                1.4)
        getattr(w,'import')(keyspdf,ROOT.RooFit.RecycleConflictNodes())

        #save fit result plot
        #showFitResult(fitVar=varName,
        #            data=fullsideBand,
        #            pdf=w.pdf('QCD_%s_%s'%(varName,ch)),
        #            categs=[''],
        #            w=w,
        #            showComponents=[],
        #            rangeX=(0,400),
        #            tagTitle='QCD_%s'%ch,
        #            outDir=opt.output)


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

        if opt.fitType<4:
            for c in ['cor','wro']:
                pdfDef=[x.format(cat) for x in  ALLPDFS[('S_%s'%c,varName)] ]
                for p in pdfDef:
                    if varName=='mjj' and opt.fitType==0: p=p.replace('(mjj','(scaledmjj')
                    pdf=w.factory(p)
                iter = pdf.getParameters(data).createIterator()
                iparam = iter.Next()
                while iparam :
                    try:
                        mcval=ws.var(iparam.GetName()).getVal()
                        iparam.setVal(mcval)
                        iparam.setConstant(True)
                    except:
                        pass
                    iparam = iter.Next()

            resFracName='f_s{1}_{0}'.format(cat,varName)
            w.factory('SUM::S_{1}_{0}( {2}[0,1]*S_cor{1}_{0}, S_wro{1}_{0} )'.format(cat,varName,resFracName))
            w.var(resFracName).setVal(1.0 if opt.fitType==3 else ws.var(resFracName).getVal() )
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
def runSimpleFit(pdf,data,poi,constr=None):

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
            print var.GetName()
            parcels.add(var)
            var = iter.Next()
        nll=ROOT.RooAddition('nllc','nllc',parcels)

    minuit=ROOT.RooMinuit(nll)
    minuit.migrad() #minimize with respect to all parameters
    minuit.minos(poi)
    r=minuit.save()
    pll=nll.createProfile(poi)
    return (r,pll)


"""
"""
def performFits(opt):

    #read workspace
    fIn=ROOT.TFile(opt.finalWS)
    w=fIn.Get('w')
    fIn.Close()

    #data to fit
    data=w.data('data')

    #fit results summary
    fitResults={}

    #run the 2D/3D fits
    if opt.fitType in [1,2]:

        for ch in ['e','mu','combined']:

            #parameter of interest
            poi=ROOT.RooArgSet()
            paramList=[] #,('eb','#varepsilon_{b}')]
            poi.add(w.var('xsec'))
            paramList.append( ('xsec','#sigma(t#bar{t})') )

            #constraints
            constr=ROOT.RooArgSet()
            constr.add( w.function('ebconstraint') )

            t='2D' if opt.fitType==1 else '3D'
            pdf=w.pdf('model_{0}_{1}'.format(ch,t))
            obsList=ROOT.RooArgSet()
            obsList.add(w.var('mjj'))
            obsList.add(w.var('mthad'))
            if t=='3D': obsList.add(w.var('mtlep'))
            #runFit(pdf,data,poi,obsList,w,opt.output)
            result,pll=runSimpleFit(pdf,data,poi,constr)
            w.saveSnapshot('fitresult_%s'%ch,pdf.getParameters(data))

            addToFitResults('%s_%s'%(ch,t),fitResults,pdf.getParameters(data))
            #pll.plotOn(xsecframe,ROOT.RooFit.Name(ch),ROOT.RooFit.ShiftToZero())

            EVENTCATEGORIES2SHOW=EVENTCATEGORIES if ch=='combined' else [ch+'1l4j'+x for x in ['2q','1b1q','2q']]
            for obs in observables:
                if obs=='mtlep' and t=='2D' : continue
                showFitResult(fitVar=obs[0],
                            data=data,
                            pdf=pdf,
                            categs=EVENTCATEGORIES2SHOW,
                            w=w,
                            showComponents=['S_cor*','S_cor*,S_wro*'],
                            rangeX=(25 if obs[0]=='mjj' else 150,
                                    300 if obs[0]=='mjj' else 400),
                            outDir=opt.output,
                            paramList=paramList,
                            pfix='_%s_%sfit'%(ch,t))

    else:
        for ch in ['e','mu','combined']:

            #parameter of interest
            poi=ROOT.RooArgSet()
            paramList=[] #,('eb','#varepsilon_{b}')]
            poi.add(w.var('xsec'))
            paramList.append( ('xsec','#sigma(t#bar{t})') )

            #constraints
            constr=ROOT.RooArgSet()
            constr.add( w.function('ebconstraint') )
            #constr.add( w.function('jsfconstraint') )

            #pdf
            pdf=w.pdf('model_%s_mjj'%ch)

            #run the fit to the mjj variable
            #obsList=ROOT.RooArgSet()
            #obsList.add(w.var('mjj'))
            #runFit(pdf,data,poi,obsList,w,opt.output)
            result,pll=runSimpleFit(pdf,data,poi,constr)
            w.saveSnapshot('fitresult_%s'%ch,pdf.getParameters(data))
            addToFitResults('%s_1D'%ch,fitResults,pdf.getParameters(data))
            #pll.plotOn(xsecframe,ROOT.RooFit.Name(ch),ROOT.RooFit.ShiftToZero())
            continue
            compsToShow=['S_cor*','S_cor*,S_wro*'] #,'S_cor*,S_wro*,W_*']
            if opt.fitType==3 :
                compsToShow=['S_cor*','S_cor*,W_*']
            if opt.fitType==4 :
                compsToShow=['S_*','S_*,W_*']
                paramList += [ ('mu_scormjj','#mu'),
                            ('sigmaL_scormjj','#sigma_{L}'),
                            ('sigmaR_scormjj','#sigma_{R}')]

            EVENTCATEGORIES2SHOW=EVENTCATEGORIES if ch=='combined' else [ch+'1l4j'+x for x in ['2q','1b1q','2q']]
            showFitResult(fitVar='mjj',
                          data=data,
                          pdf=pdf,
                          categs=EVENTCATEGORIES2SHOW,
                          w=w,
                          showComponents=compsToShow,
                          rangeX=(25,300),
                          outDir=opt.output,
                          paramList=paramList,
                          pfix='_%s_%dfit'%(ch,opt.fitType))

    printFitResults(fitResults,opt)
    w.writeToFile('finalfitworkspace_%d.root'%(opt.fitType))

    #show the likelihoods
    #xsecframe=w.var('xsec').frame(ROOT.RooFit.Bins(10),ROOT.RooFit.Range(0,100))
    #c=ROOT.TCanvas('c','c',500,500)
    #c.SetTopMargin(0.05)
    #c.SetLeftMargin(0.12)
    #c.SetRightMargin(0.03)
    #c.SetBottomMargin(0.1)
    #xsecframe.Draw()
    #xsecframe.GetXaxis().SetTitle('#sigma [nb]')
    #xsecframe.GetYaxis().SetTitle('-2 log #lambda')
    #label = ROOT.TLatex()
    #label.SetNDC()
    #label.SetTextFont(42)
    #label.SetTextSize(0.04)
    #label.DrawLatex(0.6,0.9,'#bf{CMS} #it{preliminary}')
    #c.Modified()
    #c.Update()
    #for ext in ['png','pdf']:
    #    c.SaveAs('%s/pll_%s_%dfit.%s'%(opt.output,ch,opt.fitType))
    #raw_input()

"""
"""
def addPDFToWorkspace(opt):
    #read workspace
    fIn=ROOT.TFile(opt.input)
    w=fIn.Get('w')
    fIn.Close()

    #common to all variables and/or channels
    w.factory('xsec[0,300]')
    w.factory('lumi[%f]'%lumi)
    for ch in ['e','mu']:
        w.factory('acc_%s[%f]'%(ch,acceptance[ch]))
        w.factory('eff_%s[%f]'%(ch,efficiency[ch]))
    w.factory("RooFormulaVar::ebconstraint('0.5*pow((@0-%f)/%f,2)',{eb[0.60,0.0,1.0]})"%(ebExp,ebUnc))

    w.factory("RooFormulaVar::jsfconstraint('0.5*pow(@0,2)',{jsf[-5,5]})")
    w.factory("RooFormulaVar::jerconstraint('0.5*pow(@0,2)',{jer[-5,5]})")
    w.factory("RooFormulaVar::scaledmjj('(1+@0*@1)*(1+@2*@3)*@4',{jsf,dmjj_jes,jer,dmjj_jer,mjj})")
    w.var('jsf').setVal(0.0)
    w.var('jsf').setConstant(True)
    w.var('jer').setVal(0.0)
    w.var('jer').setConstant(True)

    #instantiate PDFs
    for vname,vtit in observables:
        w.var(vname).SetTitle(vtit)

        #model QCD
        addQCDModel(w,opt,vname)

        #add signal model
        addSignalModel(w,opt,vname)

        #define the PDFs
        definePDF(w,vname)

    defineCombinedPDFs(w)
    w.writeToFile('finalworkspace.root',True)

"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True) 

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--output',    dest='output',    default='plots/Data8TeV_pp',          type='string',   help='output directory [%default]')
    parser.add_option('-i', '--input',     dest='input',     default='workspace_Data8TeV_pp.root', type='string',   help='workspace [%default]')
    parser.add_option('-s', '--signal',    dest='signal',    default='pdf_workspace_MC8.16TeV_TTbar_pPb.root', type='string',   help='signal workspace [%default]')
    parser.add_option(      '--fitType',   dest='fitType',   default=0,                  type=int,
        help='0-full signal; 1-full signal 2D; 2-full signal 3D; 3-res from MC; 4-res from CB [%default]')
    parser.add_option('-v', '--verbose',   dest='verbose',   default=0,                            type=int,        help='Verbose mode [%default]')
    parser.add_option(      '--finalWorkspace',      dest='finalWS',      default=None,            type='string',   help='final workspace to be used for the fit [%default]')

    (opt, args) = parser.parse_args()

    if opt.verbose<9 : shushRooFit()

    #create final workspace if not given
    if opt.finalWS is None: addPDFToWorkspace(opt)
    else:                   performFits(opt)


if __name__ == "__main__":
    main()
