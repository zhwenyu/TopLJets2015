#!/usr/bin/env python2.7

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
def showFitResult(ds,pdf,sampleCuts,w,showComponents=[],rangeX=(0,400)):

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
    for tag,cut in sampleCuts:
        p1.cd()
        p1.Clear()
        frame=w.var('mjj').frame(ROOT.RooFit.Range(rangeX[0],rangeX[1]))
        redData = ds.reduce(ROOT.RooFit.Cut("sample==sample::%s"%tag))
        redData.plotOn(frame)
        w.cat('sample').setLabel(tag)
        for icomp in xrange(0,len(showComponents)):
            pdf.plotOn(frame,
                       ROOT.RooFit.Components(showComponents[icomp]),
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
        for tkn,newTkn in [('eq','='),('g=','#geq'),('j','j,')]: title=title.replace(tkn,newTkn)
        label.DrawLatex(0.6,0.88,title)
        label.DrawLatex(0.6,0.84,'#chi^{2}=%3.2f'%frame.chiSquare())
        
        p2.cd()
        p2.Clear()
        hpull = frame.pullHist()
        pullFrame=w.var('mjj').frame(ROOT.RooFit.Range(rangeX[0],rangeX[1]))
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

"""
"""
def definePDF(multiplicityCategories,btagCategories,w):
    
    #signal is fit through the resonant component
    w.factory("RooFormulaVar::NSr('@0*@1',{NS[0,9999999.],fSr[0,1]})")

    #without external input assume signal is fully resonant
    w.var('fSr').setVal(1.0)
    w.var('fSr').setConstant(True)

    #global PDF
    simPDF=ROOT.RooSimultaneous('model','model',w.cat('sample'))

    #loop to define PDFs+yields per categories
    nMultCats=len(multiplicityCategories)
    fSSum,fSSumArgs='',''
    for i in xrange(0,nMultCats):

        mcName,_=multiplicityCategories[i]

        #fraction of signal in this category
        if nMultCats==1:
            w.factory('fS%s[1]'%mcName)
        else:
            if i<nMultCats-1 : 
                w.factory('fS%s[0,1]'%mcName)
                fSSum     += '@%d+'%i
                fSSumArgs += 'fS%s,'%mcName
            else : 
                w.factory("RooFormulaVar::fS%s('1.0-(%s)',{%s})"%(mcName,fSSum[:-1],fSSumArgs[:-1]))

        #b-finding efficiency
        #w.factory('eb%s[0.6,0,1]'%mcName)
        #w.factory('eb%s[1.0]'%mcName)

        #model parameters
        #w.factory('rsigma%s[5,60]'%mcName)

        #common signal parameters
        w.factory('srmu%s[80,50.,90.]'%mcName)
        w.factory('srsigma%s[10,5.,20.]'%mcName)

        #signal,background yields and total PDF
        for bcName,_ in btagCategories:
            tag=mcName+bcName
            
            #signal yields
            #expr="RooFormulaVar::NSr%s('@0*"%tag
            #if bcName=='eq0b':    expr += "TMath::Power(1-@1,2)"
            #elif  bcName=='eq1b': expr += "2*@1*(1-@1)"
            #else :                expr += "TMath::Power(@1,2)"
            #expr += "*@2',{fS%s,eb%s,NSr})"%(mcName,mcName)
            #w.factory(expr)

            #SIGNAL : fr.CB + fnov.NOV + (1-fr-fnov).EXP


            #resonant component (resolution may depend on category)
            w.factory('fr%s[0,1]'%tag)
            w.factory('sralpha%s[1.2,0.1,5.]'%tag)
            w.factory('srn%s[1.0,0.05,8.]'%tag)
            w.factory('RooCBShape::Sr%s(mjj,srmu%s,srsigma%s,sralpha%s,srn%s)'%(tag,mcName,mcName,tag,tag))

            #"turn-on" component
            w.factory('fnov%s[0,1]'%tag)
            w.factory('bnrpeak%s[30,10,50]'%tag)
            w.factory('bnrsigma%s[10,50]'%tag)
            w.factory('bnrtail%s[-10,0]'%tag)
            w.factory("RooNovosibirsk::Bnov%s(mjj,bnrpeak%s,bnrsigma%s,bnrtail%s)"%(tag,tag,tag,tag))
            
            #exponential component
            #w.factory("RooExponential::Bexp%s(mjj,bnrlambda%s[-0.5,-10,0])"%(tag,tag))
            
            #category PDF :
            #pdf=w.factory('SUM:model%s(fr%s*Sr%s,fnov%s*Bnov%s,Bexp%s)'%(tag,tag,tag,tag,tag,tag))
            pdf=w.factory('SUM:model%s(fr%s*Sr%s,Bnov%s)'%(tag,tag,tag,tag))

            #add to simultaneous pdf
            simPDF.addPdf(w.pdf('model%s'%tag),tag)

    #import simultaneous pdf to workspace
    getattr(w,'import')(simPDF,ROOT.RooFit.RecycleConflictNodes())
    return w.pdf(simPDF.GetName())


"""
defines the categories to use
"""
def defineCategories(categorizationMode,dataName,w):

    multiplicityCategories,btagCategories=None,None
    sampleCuts=[]
    if categorizationMode==1:
        #multiplicityCategories = [('eq3j','nj==3'),('geq4j','nj>=4')]
        multiplicityCategories = [('geq4j','nj>=4')]
        btagCategories         = [('eq0b','nb==0'),('eq1b','nb==1'),('geq2b','nb>=2')]
    

    #add categories to the workspace
    sample=ROOT.RooCategory('sample','sample')
    sample.defineType('fail')
    for mcName,mcCut in multiplicityCategories:
        for bcName,bcCut in btagCategories:
            tag=mcName+bcName
            sample.defineType(tag)
            sampleCuts.append( (tag, mcCut + ' && ' + bcCut) )
    getattr(w,'import')(sample)

    #all done
    return multiplicityCategories,btagCategories,sampleCuts

"""
categorizes the dataset according to the option chosen
"""
def categorizeDataSet(sampleCuts,dataName,w,bkgSF=0.,lumiSF=100.):
    
    data=w.data(dataName)

    #extend argumments of the dataset with a category
    args=ROOT.RooArgSet()
    args.add( w.cat('sample') )   
    args.add( w.var('w') )  
    entryIt=data.get().createIterator()
    var = entryIt.Next()
    while var :
        args.add( w.var( var.GetName() ) )
        var = entryIt.Next()

    #start dataset
    combData=ROOT.RooDataSet('comb%s'%dataName,'comb%s'%dataName, args, 'w')

    #fill dataset by slicing the original one
    for tag,cut in sampleCuts:
        redData=data.reduce(cut + ' && id>=0')
        for i in xrange(0,redData.numEntries()):
            entry=redData.get(i)
            args.setCatLabel('sample',tag)
           
            entryIt = entry.createIterator()
            var = entryIt.Next()
            while var :
                args.find( var.GetName() ).setVal( var.getVal() )
                var = entryIt.Next()

            idVal=args.find('id').getVal()
            wgt=redData.weight()*lumiSF
            if idVal>0 : wgt*=bkgSF
            args.find('w').setVal( wgt )
            
            combData.add( args, wgt )

    getattr(w,'import')(combData)
    return w.data('comb%s'%dataName)

"""
minimizes the likelihood and profiles the poi
"""
def runFit(pdf,data,poi,xmin,xmax):
    nll=pdf.createNLL(data,ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Range(xmin,xmax),ROOT.RooFit.NumCPU(2))
    ROOT.RooMinimizer(nll).migrad() #minimize with respect to all parameters    
    pll=nll.createProfile(poi)
    return (nll,pll)


"""
"""
def runPseudoExperiments(sampleCuts,origData,poi,rangeX,w,nPexp):

        #arguments for the combined dataset
        args=ROOT.RooArgSet()
        args.add( w.cat('sample') ) 
        args.add( w.var('mjj') )

        #observable to fit
        obsArgs=ROOT.RooArgSet()
        obsArgs.add( w.var('mjj') )
                            
        for i in xrange(0,nPexp):
                
            combData=ROOT.RooDataSet('pe','pe', args)            

            #generate the dataset for each category from the best fit to MC
            w.loadSnapshot('bestfit')
            for tag,cut in sampleCuts:
                args.setCatLabel('sample',tag)
                pdf=w.pdf('model'+tag)

                #generate randomly from a Poisson centred at the expected # of events
                nEvts=ROOT.gRandom.Poisson( origData.sumEntries(cut) )
                ds=pdf.generate(obsArgs,nEvts)
            
                #add events to the combined dataset
                for j in xrange(0,ds.numEntries()):
                    entry=ds.get(j)
                    args.find('mjj').setVal( entry.find('mjj').getVal() )
                    combData.add(args)

                #all done with this one
                ds.Delete()

            nll,pll=runFit(pdf=w.pdf('model'), data=combData,poi=poi, xmin=rangeX[0], xmax=rangeX[1])

            #for debug purposes only
            #showFitResult(ds=combData,pdf=w.pdf('model'),sampleCuts=sampleCuts,w=w,rangeX=(30,200))

            #free mem
            combData.Delete()
            nll.Delete()
            pll.Delete()
            

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',     dest='input',     help='input file [%default]',              default='dijetWorkspace.root',  type='string')
    parser.add_option(      '--data',      dest='data',      help='dataset name [%default]',            default='mc',                   type='string')
    parser.add_option(      '--catMode',   dest='catMode',   help='categorization mode [%default]',     default=1,                      type='int')
    parser.add_option(      '--bkgSF',     dest='bkgSF',     help='background scale factor [%default]', default=1.0,                    type='float')
    parser.add_option(      '--lumiSF',    dest='lumiSF',    help='lumi scale factor [%default]',       default=1.0,                    type='float')
    parser.add_option(      '--fitRange',  dest='fitRange',  help='fit range [%default]',               default='30,200',               type='string')
    parser.add_option(      '--runPE',     dest='runPE',     help='#p.e. to run [%default]',            default=0,                      type='int')
    parser.add_option('-v', '--verbose',   dest='verbose',   default=0, type=int,   help='Verbose mode [%default]')
    (opt, args) = parser.parse_args()

    if opt.verbose<9 : shushRooFit()
     
    #read workspace
    fIn=ROOT.TFile(opt.input)
    w=fIn.Get('w')

    #define the categories
    multiplicityCategories,btagCategories,sampleCuts=defineCategories(opt.catMode,opt.data,w)
    print('Defined %d exclusive categories'%len(sampleCuts))

    #get the categorized dataset
    combData=categorizeDataSet(sampleCuts=sampleCuts,dataName=opt.data,w=w,bkgSF=opt.bkgSF,lumiSF=opt.lumiSF)
    print('Combined dataset to fit has %d entries'%combData.sumEntries())

    #define the PDFs
    combPDF=definePDF(multiplicityCategories,btagCategories,w)
    print('Defined combined PDF')

    #define all variables and the POI
    allVars=w.allVars()
    poi=ROOT.RooArgSet()
    varIt=allVars.createIterator()
    var = varIt.Next()
    while var :
        varName=var.GetName()
        if varName.find('srmu')==0 or varName.find('fr')==0 :
            poi.add( var )
        var=varIt.Next()

    #perform simultaneous fit - profile all parameters except the POI
    xmin,xmax=[float(x) for x in opt.fitRange.split(',')]
    runFit(pdf=combPDF, data=combData,poi=poi, xmin=xmin, xmax=xmax)

    #save values postfit
    w.saveSnapshot("poi",poi)
    w.saveSnapshot("bestfit",allVars)
    
    #show results
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    print 'Saving fit result plots'
    showFitResult(ds=combData,pdf=combPDF,sampleCuts=sampleCuts,w=w,rangeX=(xmin,xmax),showComponents=['B*'])

    if opt.runPE>0 :
        print 'Will run %d pseudo-experiments'%opt.runPE
        runPseudoExperiments(sampleCuts=sampleCuts,origData=combData,poi=poi,rangeX=(xmin,xmax),w=w,nPexp=opt.runPE)



if __name__ == "__main__":
    main()

