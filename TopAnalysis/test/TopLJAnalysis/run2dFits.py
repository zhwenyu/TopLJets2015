#!/usr/bin/env python2.7

#from prepareWorkspace import EVENTCATEGORIES
EVENTCATEGORIES=[
    '1f4j1q', '1l4j2b','1l4j1b','1l4j1q',
    #'1l3j1b','1l3j1q','1f3j1q',
    ]

FITVARS=['m_wjj','m_bwjj']
XRANGES={'m_wjj':(0,400),'m_bwjj':(60,400)}

import ROOT
import optparse
import os,sys

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
display the results of the fit
"""
def showFitResult(fitVar,data,pdf,categs,w,showComponents=[],rangeX=(0,400),postfix='',outDir='plots/'):

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0)
    c.SetLeftMargin(0)
    c.SetRightMargin(0)
    c.SetBottomMargin(0)
    p1 = ROOT.TPad('p1','p1',0.0,0.85,1.0,0.0)
    p1.SetRightMargin(0.05)
    p1.SetLeftMargin(0.12)
    p1.SetTopMargin(0.008)
    p1.SetBottomMargin(0.12)
    p1.Draw()
    c.cd()
    p2 = ROOT.TPad('p2','p2',0.0,0.85,1.0,1.0)
    p2.SetBottomMargin(0.005)
    p2.SetRightMargin(0.05)
    p2.SetLeftMargin(0.12)
    p2.SetTopMargin(0.05)
    p2.SetGridy(True)
    p2.Draw()

    #show fit result
    for tag in categs:
        p1.cd()
        p1.Clear()
        frame=w.var(fitVar).frame(ROOT.RooFit.Range(rangeX[0],rangeX[1]),ROOT.RooFit.Bins(50))
        redData = data.reduce(ROOT.RooFit.Cut("sample==sample::%s"%tag))
        redData.plotOn(frame)
        w.cat('sample').setLabel(tag)
        if pdf:
            for icomp in xrange(0,len(showComponents)):
                comps=showComponents[icomp]
                if comps=='S_*':
                    color=ROOT.TColor.GetColor('#d73027')
                    pdf.plotOn(frame,
                            ROOT.RooFit.Components(comps),
                            ROOT.RooFit.Slice(w.cat('sample'),tag),
                            ROOT.RooFit.ProjWData(redData),
                            ROOT.RooFit.LineColor(color),
                            ROOT.RooFit.FillColor(color),
                            ROOT.RooFit.FillStyle(1001),
                            ROOT.RooFit.DrawOption('f'),
                            ROOT.RooFit.LineWidth(1))
                else:
                     pdf.plotOn(frame,
                            ROOT.RooFit.Components(comps),
                            ROOT.RooFit.Slice(w.cat('sample'),tag),
                            ROOT.RooFit.ProjWData(redData),
                            ROOT.RooFit.LineStyle(2),
                            ROOT.RooFit.LineColor(ROOT.kGray+icomp+1),
                            ROOT.RooFit.LineWidth(1))
            pdf.plotOn(frame,
                    ROOT.RooFit.Slice(w.cat('sample'),tag),
                    ROOT.RooFit.ProjWData(redData),
                    ROOT.RooFit.LineColor(ROOT.kBlue),
                    ROOT.RooFit.LineWidth(2))
        frame.Draw()
        frame.GetYaxis().SetTitle("Entries")
        frame.GetYaxis().SetTitleOffset(1.0)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetLabelSize(0.04)
        xtitle="M(jj) [GeV]" if fitVar=='m_wjj' else 'M(bjj) [GeV]'
        frame.GetXaxis().SetTitle(xtitle)

        label = ROOT.TLatex()
        label.SetNDC()
        label.SetTextFont(42)
        label.SetTextSize(0.04)
        label.DrawLatex(0.6,0.92,'#bf{CMS} #it{preliminary}')
        tagTitle ='iso' if '1l' in tag else 'non-iso'
        tagTitle+=',#geq4j' if '4j' in tag else ',=3j'
        tagTitle+=','+tag[-2:].replace('q','non-b')
        label.DrawLatex(0.6,0.88,tagTitle)
        label.DrawLatex(0.6,0.84,'#chi^{2}=%3.2f'%frame.chiSquare())
        j='3j' if '3j' in tag else '4j'
        ivar=0
        for var,tit in [('Nsig','N_{cp}(t#bar{t})'),
                        ('sig_mu_%s'%fitVar,'#mu'),
                        ('sig_sigma_%s'%fitVar,'#sigma'),
                        ('eb_%s'%j,'#epsilon_{b}')]:
            label.DrawLatex(0.6,
                            0.80-0.04*ivar,
                            '#scale[0.8]{%s=%3.2f#pm%3.2f}'%(tit,
                                                             w.var(var).getVal(),
                                                             w.var(var).getError()) )
            ivar+=1

        p2.cd()
        p2.Clear()
        if pdf:
            hpull = frame.pullHist()
            pullFrame=w.var(fitVar).frame(ROOT.RooFit.Range(rangeX[0],rangeX[1]))
            pullFrame.addPlotable(hpull,"P") ;
            pullFrame.Draw()
            pullFrame.GetYaxis().SetTitle("Pull")
            pullFrame.GetYaxis().SetTitleSize(0.25)
            pullFrame.GetYaxis().SetLabelSize(0.25)
            pullFrame.GetXaxis().SetTitleSize(0)
            pullFrame.GetXaxis().SetLabelSize(0)
            pullFrame.GetYaxis().SetTitleOffset(0.15)
            pullFrame.GetYaxis().SetNdivisions(4)
            pullFrame.GetYaxis().SetRangeUser(-3.1,3.1)
            pullFrame.GetXaxis().SetTitleOffset(0.8)

        c.cd()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('%s/%s_%s_%sfit.%s'%(outDir,fitVar,tag,postfix,ext))

"""
"""
def definePDF(w,fitVars=['m_wjj','m_wbjj'],paramQCD=False,procList=['S','W','QCD']):



    #EFFICIENCIES
    #efficiencies will control the yields accross categories through binominal expansion
    w.factory('e4jtoe3j[1.0]')
    for j in ['3j','4j'] :
        w.factory('eb_%s[0.8,0.3,0.9]'%j)
        if not 'S' in procList:
            w.var('eb_%s'%j).setVal(0)
            w.var('eb_%s'%j).setConstant(True)

    #YIELDS
    #signal yields are parametrized as function of the b-finding efficiencies
    #background yields are category-dependent like
    w.factory("Nsig[1000,0,100000]")
    if not 'S' in procList:
        w.var('Nsig').setVal(0)
        w.var('Nsig').setConstant(True)
    for cat in EVENTCATEGORIES:
        j='3j' if '3j' in cat  else '4j'
        expr='@0'
        if j=='3j':
            if '1q' in cat : expr += '*(1-@1)'
            else           : expr += '*@1'
            expr += '*(1-@2)'
        else:
            if '1q' in cat   : expr += '*pow(1-@1,2)'
            elif '1b' in cat : expr += '*2*@1*(1-@1)'
            else             : expr += '*pow(@1,2)'
            expr += '*@2'
        if '1f' in cat:
            w.factory("RooFormulaVar::Nsig_%s('0*@0',{Nsig})"%(cat))
        else:
            w.factory("RooFormulaVar::Nsig_%s('%s',{Nsig,eb_%s,e4jtoe3j})"%(cat,expr,j))

        #ignore w-like in a pure QCD region
        w.factory('Nbkg_%s[1000,0,100000]'%cat)
        if '1f' in cat and '1q' in cat:
            w.factory('fqcd_%s[1.0]'%cat)
        else:
            if '2b'   in cat : w.factory('fqcd_%s[0.0,0.2]'%cat)
            elif '1b' in cat : w.factory('fqcd_%s[0.0,0.0,0.5]'%cat)
            else             : w.factory('fqcd_%s[0.0,0.0,0.5]'%cat)
            if not 'W' in procList:
                w.factory('fqcd_%s'%cat).setVal(1)
                w.factory('fqcd_%s'%cat).setConstant(True)
        w.factory("RooFormulaVar::Nqcd_%s('@0*@1',{Nbkg_%s,fqcd_%s})"%(cat,cat,cat))
        w.factory("RooFormulaVar::Nw_%s('@0*(1-@1)',{Nbkg_%s,fqcd_%s})"%(cat,cat,cat))

    #
    # PDFs per variable
    #
    for fitVar in fitVars:

        #global PDF
        fitVarSimPDF=ROOT.RooSimultaneous('%s_pdf'%fitVar,'%s_pdf'%fitVar,w.cat('sample'))

        #signal model is a crystal ball with common mean and width accross categories
        if fitVar=='m_wjj':
            w.factory('sig_mu_%s[80,0,150]'%fitVar)
            w.factory('sig_sigma_%s[1,20]'%fitVar)
        else:
            w.factory('sig_mu_%s[175,150,200]'%fitVar)
            w.factory('sig_sigma_%s[1,20]'%fitVar)
        for j in ['3j','4j'] :
            w.factory('sig_alpha_%s_%s[0.01,5]'%(fitVar,j))
            w.factory('sig_n_%s_%s[0.05,3]'%(fitVar,j))
            w.factory('RooCBShape::S_%s_%s(%s,sig_mu_%s,sig_sigma_%s,sig_alpha_%s_%s,sig_n_%s_%s)'%(fitVar,j,fitVar,fitVar,fitVar,fitVar,j,fitVar,j))

        #QCD background is only jet dependent from RooKeysPDF if non-parametric, otherwise use simple functions
        for j in ['3j','4j'] :
            if not paramQCD:
                qcdCat='1f%s1q'%j
                redData = w.data('data').reduce(ROOT.RooFit.Cut("sample==sample::%s"%qcdCat))
                #if redData.numEntries()>500:
                #    print 'reducing numentries to 500'
                #    redData=redData.reduce(ROOT.RooFit.EventRange(0,750))
                keyspdf=ROOT.RooKeysPdf("QCD_%s_%s"%(fitVar,j),
                                        "QCD_%s_%s"%(fitVar,j),
                                        w.var(fitVar),
                                        redData,ROOT.RooKeysPdf.NoMirror,1) #Both)
                getattr(w,'import')(keyspdf,ROOT.RooFit.RecycleConflictNodes())
            else:
                if fitVar=='m_wjj':
                    w.factory('qcd_alpha_%s_%s[10,0,100]'%(fitVar,j))
                    w.factory('qcd_beta_%s_%s[0,0,1]'%(fitVar,j))
                    w.factory("EXPR::QCD_%s_%s('@0*TMath::Exp(-0.5*TMath::Power(@0/(@1+@2*@0),2))',{%s,qcd_alpha_%s_%s,qcd_beta_%s_%s})"%(fitVar,j,fitVar,fitVar,j,fitVar,j))
                    #w.factory('qcd_peak_%s_%s[0,50]'%(fitVar,j))
                    #w.factory('qcd_sigma_%s_%s[10,50]'%(fitVar,j))
                    #w.factory("RooLandau::QCD_%s_%s(%s,qcd_peak_%s_%s,qcd_sigma_%s_%s)"%(fitVar,j,fitVar,fitVar,j,fitVar,j))
                    #w.factory('qcd_tail_%s_%s[-10,0]'%(fitVar,j))
                    #w.factory("RooNovosibirsk::QCD_%s_%s(%s,qcd_peak_%s_%s,qcd_sigma_%s_%s,qcd_tail_%s_%s)"%(fitVar,j,fitVar,fitVar,j,fitVar,j,fitVar,j))
                else:
                    w.factory('qcd_alpha_%s_%s[10,0,100]'%(fitVar,j))
                    w.factory('qcd_beta_%s_%s[0,0,1]'%(fitVar,j))
                    w.factory("EXPR::QCD_%s_%s('@0*TMath::Exp(-0.5*TMath::Power(@0/(@1+@2*@0),2))',{%s,qcd_alpha_%s_%s,qcd_beta_%s_%s})"%(fitVar,j,fitVar,fitVar,j,fitVar,j))
                    #w.factory('qcd_peak_%s_%s[20,200]'%(fitVar,j))
                    #w.factory('qcd_sigma_%s_%s[10,200]'%(fitVar,j))
                    #w.factory("RooLandau::QCD_%s_%s(%s,qcd_peak_%s_%s,qcd_sigma_%s_%s)"%(fitVar,j,fitVar,fitVar,j,fitVar,j))

        #W-like is category dependent
        for cat in EVENTCATEGORIES:
            if fitVar=='m_wjj':
                w.factory('w_alpha_%s_%s[10,0,100]'%(fitVar,cat))
                w.factory('w_beta_%s_%s[0,0,1]'%(fitVar,cat))
                w.factory("EXPR::W_%s_%s('@0*TMath::Exp(-0.5*TMath::Power(@0/(@1+@2*@0),2))',{%s,w_alpha_%s_%s,w_beta_%s_%s})"%(fitVar,cat,fitVar,fitVar,cat,fitVar,cat))
            else:
                w.factory('w_peak_%s_%s[10,200]'%(fitVar,cat))
                w.factory('w_sigma_%s_%s[10,200]'%(fitVar,cat))
                w.factory("RooLandau::W_%s_%s(%s,w_peak_%s_%s,w_sigma_%s_%s)"%(fitVar,cat,fitVar,fitVar,cat,fitVar,cat))

        #combined model
        for cat in EVENTCATEGORIES:
            j='3j' if '3j' in cat else '4j'
            w.factory('SUM:model_%s_%s(Nsig_%s*S_%s_%s,Nqcd_%s*QCD_%s_%s,Nw_%s*W_%s_%s)'%(fitVar,cat,cat,fitVar,j,cat,fitVar,j,cat,fitVar,cat))
            fitVarSimPDF.addPdf(w.pdf('model_%s_%s'%(fitVar,cat)),cat)

        #import simultaneous pdf to workspace
        getattr(w,'import')(fitVarSimPDF,ROOT.RooFit.RecycleConflictNodes())


    #now define the 2D pdf mjj vs mbjj for each process
    simPDF=ROOT.RooSimultaneous('pdf','pdf',w.cat('sample'))
    for proc in ['S','QCD','W']:
        categs = EVENTCATEGORIES if proc=='W' else ['3j','4j']
        for cat in categs:
            expr='PROD:%s_%s('%(proc,cat)
            for fitVar in fitVars:
                expr += '%s_%s_%s,'%(proc,fitVar,cat)
            expr=expr[:-1]+')'
            w.factory(expr)

    #sum up the 2D pdf in the categories
    for cat in EVENTCATEGORIES:
        j='3j' if '3j' in cat else '4j'
        pdf=w.factory('SUM:model_%s(Nsig_%s*S_%s,Nqcd_%s*QCD_%s,Nw_%s*W_%s)'%(cat,cat,j,cat,j,cat,cat))
        simPDF.addPdf(w.pdf('model_%s'%cat),cat)
    getattr(w,'import')(simPDF,ROOT.RooFit.RecycleConflictNodes())


"""
minimizes the likelihood and profiles the poi
"""
def runFit(pdf,data,poi):
    nll=pdf.createNLL(data,
                      ROOT.RooFit.Extended(True),
                      ROOT.RooFit.SumW2Error(True),
                      ROOT.RooFit.NumCPU(2))
    ROOT.RooMinimizer(nll).migrad() #minimize with respect to all parameters
    pll=nll.createProfile(poi)
    return (nll,pll)

"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(False)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--output',    dest='output',    default='plots/Data8TeV_pp',          type='string',   help='output directory [%default]')
    parser.add_option('--paramQCD', dest='paramQCD', default=False, action='store_true', help='do not use RooKeysPdf for QCD [%default]')
    parser.add_option('--procList', dest='procList', default='S,W,QCD',  help='parameterize only process in this list [%default]')
    parser.add_option('-i', '--input',     dest='input',     default='workspace_Data8TeV_pp.root', type='string',   help='workspace [%default]')
    parser.add_option('-v', '--verbose',   dest='verbose',   default=0,                            type=int,        help='Verbose mode [%default]')
    (opt, args) = parser.parse_args()

    if opt.verbose<9 : shushRooFit()

    #read workspace
    fIn=ROOT.TFile(opt.input)
    w=fIn.Get('w')

    #data to fit
    data=w.data('data')

    #define the PDFs
    definePDF(w,fitVars=FITVARS,paramQCD=opt.paramQCD,procList=opt.procList.split(','))

    #perform simultaneous fit - profile all parameters except the POI
    #for fitVar in fitVars:
    #    xmin,xmax=xRanges[fitVar]
    #    pdf=w.pdf('%s_pdf'%fitVar)

    #    poi=ROOT.RooArgSet()
    #    poiVar='sig_mu_%s'%fitVar
    #    poi.add(w.var(poiVar))
    #    ll,pll=runFit(pdf=pdf, data=data,poi=poi)
    #    print 'Likelihood minimized for',fitVar

    #    showFitResult(fitVar=fitVar,data=data,pdf=pdf,categs=EVENTCATEGORIES,w=w,showComponents=['S_*','S_*,W_*'],rangeX=(xmin,xmax),outDir=opt.output)
    #    print 'Fit result plots saved'

    #final 2d fit
    pdf=w.pdf('pdf')
    poi=ROOT.RooArgSet()
    poi.add(w.var('Nsig'))
    ll,pll=runFit(pdf=pdf, data=data,poi=poi)
    w.saveSnapshot('fit',w.allVars())
    for fitVar in FITVARS:
        xmin,xmax=XRANGES[fitVar]
        showFitResult(fitVar=fitVar,data=data,pdf=pdf,categs=EVENTCATEGORIES,w=w,showComponents=['S_*','S_*,W_*'],rangeX=(xmin,xmax),postfix='2d',outDir=opt.output)
    print 'Fit result plots saved'
    w.writeToFile(opt.input.replace('.root','_fit.root'))


if __name__ == "__main__":
    main()
