#!/usr/bin/env python2.7

#from prepareWorkspace import EVENTCATEGORIES
EVENTCATEGORIES=[
    '1l4j2b','1l4j1b','1l4j1q','1f4j1q','1f4j1b'#,'1f4j2b',
    #'1l3j1b','1l3j1q','1f3j1q',#'1f3j1b'
    ]
    
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
def showFitResult(fitVar,data,pdf,categs,w,showComponents=[],rangeX=(0,400)):

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
        frame.GetXaxis().SetTitle("M(j,j) [GeV]")
        label = ROOT.TLatex()
        label.SetNDC()
        label.SetTextFont(42)
        label.SetTextSize(0.04)
        label.DrawLatex(0.6,0.92,'#bf{CMS} #it{preliminary}')
        title=tag
        #for tkn,newTkn in [('eq','='),('g=','#geq'),('j','j,')]: title=title.replace(tkn,newTkn)
        label.DrawLatex(0.6,0.88,title)
        #label.DrawLatex(0.6,0.84,'#chi^{2}=%3.2f'%frame.chiSquare())
        
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
        c.SaveAs('%s_fit.png'%tag)
        raw_input()
        
"""
"""
def definePDF(w,fitVar):

    #global PDF
    simPDF=ROOT.RooSimultaneous('pdf','pdf',w.cat('sample'))
    
    #signal model is a crystal ball with common mean and width accross categories
    #background models are Novosibirks functions
    if fitVar=='m_wjj':
        w.factory('sig_mu[0,300]')
        w.factory('sig_sigma[5,40]')
    else:
        w.factory('sig_mu[160,180]')
        w.factory('sig_sigma[10,40]')
    w.factory('e4jtoe3j[0.5,0,1]')
    w.factory('el[0.9,0.4,1]')
    for j in ['3j','4j'] :
        w.factory('eb_%s[0.5,0.4,1]'%j)
        w.factory('sig_alpha_%s[1.2,0.1,5.]'%j)
        w.factory('sig_n_%s[1.0,0.05,8.]'%j)
        w.factory('RooCBShape::S_%s(%s,sig_mu,sig_sigma,sig_alpha_%s,sig_n_%s)'%(j,fitVar,j,j))

        if fitVar=='m_wjj':
            w.factory('qcd_peak%s[30,10,50]'%j)
            w.factory('qcd_sigma%s[10,50]'%j)
        else:
            w.factory('qcd_peak%s[100,50,150]'%j)
            w.factory('qcd_sigma%s[10,80]'%j)
        w.factory('qcd_tail%s[-10,0]'%j)
        w.factory("RooNovosibirsk::QCD_%s(%s,qcd_peak%s,qcd_sigma%s,qcd_tail%s)"%(j,fitVar,j,j,j))

        if fitVar=='m_wjj':
            w.factory('w_peak%s[30,10,50]'%j)
            w.factory('w_sigma%s[10,50]'%j)
        else:
            w.factory('w_peak%s[100,50,150]'%j)
            w.factory('w_sigma%s[10,80]'%j)
        w.factory('w_tail%s[-10,0]'%j)
        w.factory("RooNovosibirsk::W_%s(%s,w_peak%s,w_sigma%s,w_tail%s)"%(j,fitVar,j,j,j))
        
    #signal yields are parametrized as lepton isolation and b-finding efficiencies
    #background yields are category-dependent like
    w.factory("Nsig[0,1000]")
    for cat in EVENTCATEGORIES:
        j='3j' if '3j' in cat  else '4j'
        expr='@0'
        if '1f' in cat : expr += '*(1-@1)'
        else           : expr += '*@1'
        if j=='3j':
            if '1q' in cat : expr += '*(1-@2)'
            else           : expr += '*@2'
            expr += '*(1-@3)'
        else:
            if '1q' in cat   : expr += '*pow(1-@2,2)'
            elif '1b' in cat : expr += '*2*@2*(1-@2)'
            else             : expr += '*pow(@2,2)'
            expr += '*@3'
        w.factory("RooFormulaVar::Nsig_%s('%s',{Nsig,el,eb_%s,e4jtoe3j})"%(cat,expr,j))
        
        w.factory('Nqcd_%s[0,10000]'%cat)
        if '1f' in cat and '1q' in cat :
            w.factory('Nw_%s[0]'%cat)
        else:
            w.factory('Nw_%s[0,10000]'%cat)
        
        pdf=w.factory('SUM:model_%s(Nsig_%s*S_%s,Nqcd_%s*QCD_%s,Nw_%s*W_%s)'%(cat,cat,j,cat,j,cat,j))
        simPDF.addPdf(w.pdf('model_%s'%cat),cat)

    #import simultaneous pdf to workspace
    getattr(w,'import')(simPDF,ROOT.RooFit.RecycleConflictNodes())
    return w.pdf(simPDF.GetName())

"""
"""
def defineTemplatePDF(dataMC,w,bins,xmin,xmax,fitVar):

    #global PDF
    simPDF=ROOT.RooSimultaneous('pdf','pdf',w.cat('sample'))
    
    #signal model is a crystal ball with common mean and width accross categories
    #background models are Novosibirks functions
    if fitVar=='m_wjj':
        w.factory('sig_mu[-0.0,-2.0,2.0]')
        w.factory('sig_sigma[1,0.7,10.0]')
    else:
        w.factory('sig_mu[172.5,120,220]')
        w.factory('sig_sigma[20,5,50]')
    w.factory('e4jtoe3j[0.5,0,1]')
    w.factory('el[0.9,0.4,1]')
    for j in ['3j','4j'] :
        w.factory('eb_%s[0.5,0.4,1]'%j)
        w.factory('sig_alpha_%s[1.2,0.1,5.]'%j)
        w.factory('sig_n_%s[1.0,0.05,8.]'%j)

        if fitVar=='m_wjj':
            w.factory('qcd_peak%s[30,10,50]'%j)
            w.factory('qcd_sigma%s[10,50]'%j)
            w.factory('qcd_tail%s[-10,0]'%j)
            w.factory("RooNovosibirsk::QCD_%s(%s,qcd_peak%s,qcd_sigma%s,qcd_tail%s)"%(j,fitVar,j,j,j))
        else:
            w.factory('qcd_peak%s[140,20,200]'%j)
            w.factory('qcd_sigma%s[50,10,200]'%j)
            w.factory("RooLandau::QCD_%s(%s,qcd_peak%s,qcd_sigma%s)"%(j,fitVar,j,j))

        if fitVar=='m_wjj':
            w.factory('w_peak%s[30,10,50]'%j)
            w.factory('w_sigma%s[10,50]'%j)
            w.factory('w_tail%s[-10,0]'%j)
            w.factory("RooNovosibirsk::W_%s(%s,w_peak%s,w_sigma%s,w_tail%s)"%(j,fitVar,j,j,j))
        else:
            w.factory('w_peak%s[140,10,200]'%j)
            w.factory('w_sigma%s[50,10,200]'%j)
            w.factory("RooLandau::W_%s(%s,w_peak%s,w_sigma%s)"%(j,fitVar,j,j))
        
    #signal yields are parametrized as lepton isolation and b-finding efficiencies
    #background yields are category-dependent like
    w.factory("Nsig[0,1000]")
    for cat in EVENTCATEGORIES:
        j='3j' if '3j' in cat  else '4j'
        expr='@0'
        if '1f' in cat : expr += '*(1-@1)'
        else           : expr += '*@1'
        if j=='3j':
            if '1q' in cat : expr += '*(1-@2)'
            else           : expr += '*@2'
            expr += '*(1-@3)'
        else:
            if '1q' in cat   : expr += '*pow(1-@2,2)'
            elif '1b' in cat : expr += '*2*@2*(1-@2)'
            else             : expr += '*pow(@2,2)'
            expr += '*@3'
        w.factory("RooFormulaVar::Nsig_%s('%s',{Nsig,el,eb_%s,e4jtoe3j})"%(cat,expr,j))
        
        w.factory('Nqcd_%s[0,10000]'%cat)
        if '1f' in cat and '1q' in cat :
            w.factory('Nw_%s[0]'%cat)
        else:
            w.factory('Nw_%s[0,10000]'%cat)
        
        redDataMC = ROOT.RooDataSet(dataMC.reduce(ROOT.RooFit.Cut("sample==sample::%s"%cat)))
        redDataMC.SetTitle("%s_%s"%(dataMC.GetTitle(),cat))
        redDataMC.SetName("%s_%s"%(dataMC.GetName(),cat))
        
        var = ROOT.RooRealVar(fitVar,fitVar, xmin,xmax);
        var.setBins(bins)
        template = ROOT.RooDataHist(redDataMC.GetName(),redDataMC.GetTitle(),ROOT.RooArgSet(var),redDataMC)
        getattr(w,'import')(template)
        
                
        w.factory("HistPdf::Sig_%s(%s,%s)"%(cat,fitVar,template.GetName()))
        w.factory("Gaussian::SigRes(%s,sig_mu,sig_sigma)"%fitVar)
        w.factory("FCONV::S_%s(%s, Sig_%s, SigRes)"%(cat,fitVar,cat));
        
        pdf=w.factory('SUM:model_%s(Nsig_%s*S_%s,Nqcd_%s*QCD_%s,Nw_%s*W_%s)'%(cat,cat,cat,cat,j,cat,j))
        simPDF.addPdf(w.pdf('model_%s'%cat),cat)

    #import simultaneous pdf to workspace
    getattr(w,'import')(simPDF,ROOT.RooFit.RecycleConflictNodes())
    return w.pdf(simPDF.GetName())



"""
minimizes the likelihood and profiles the poi
"""
def runFit(pdf,data,poi,xmin,xmax):
    nll=pdf.createNLL(data,
                      ROOT.RooFit.Extended(True),
                      ROOT.RooFit.SumW2Error(True),
                      ROOT.RooFit.Range(xmin,xmax),
                      ROOT.RooFit.NumCPU(2))
    ROOT.RooMinimizer(nll).migrad() #minimize with respect to all parameters    
    pll=nll.createProfile(poi)
    return (nll,pll)

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-v', '--verbose',   dest='verbose',   default=0, type=int,   help='Verbose mode [%default]')
    parser.add_option('-f','--fit',   dest='fit',   default=True,   help='Enable template(binned) fit [%default]')
    (opt, args) = parser.parse_args()

    if opt.verbose<9 : shushRooFit()

    #hard-coded : move to option parser
    fitVar='m_wjj'
    #fitVar='m_bwjj'
    xmin,xmax=0,500

    #read workspace
    fIn=ROOT.TFile('workspace.root')
    w=fIn.Get('w')

    #data to fit
    data=w.data('data')

    if opt.fit:
        #read MC workspace
        fInMC=ROOT.TFile('workspace_MC.root')
        wMC=fInMC.Get('w')

        dataMC=wMC.data('data_withWeights')

        #define the PDFs
        bins=50
        pdf=defineTemplatePDF(dataMC,w,bins,xmin,xmax,fitVar=fitVar)
    else:
        #define the PDFs                                                                                                                                        
        pdf=definePDF(w,fitVar=fitVar)
    print 'Defined combined PDF'

    #perform simultaneous fit - profile all parameters except the POI
    poi=ROOT.RooArgSet()
    poiVar='sig_mu'
    poi.add(w.var(poiVar))
    ll,pll=runFit(pdf=pdf, data=data,poi=poi, xmin=xmin, xmax=xmax)
    print 'Likelihood minimized'

    #frame=w.var(poiVar).frame()
    #pll.plotOn(frame)
    #frame.Draw()
    #raw_input()
    
    #save values postfit
    #w.saveSnapshot("poi",poi)
    #w.saveSnapshot("bestfit",allVars)
    
    #show results
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(False)
    #print 'Saving fit result plots'
    showFitResult(fitVar=fitVar,data=data,pdf=pdf,categs=EVENTCATEGORIES,w=w,showComponents=['S_*','S_*,W_*'],rangeX=(xmin,xmax))




if __name__ == "__main__":
    main()

