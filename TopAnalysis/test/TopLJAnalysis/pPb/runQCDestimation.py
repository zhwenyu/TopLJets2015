import ROOT
import sys
import pickle

QCDNORM={}

"""
"""
def getHistos(args,var):

    #get templates from file
    histos={'qcd':{}}
    for arg in args:
        proc,url=arg.split(':')
        fIn=ROOT.TFile.Open(url)
        histos[proc]={}
        for ch in ['e','mu']:
            for cat in ['1l4j2q','1l4j1b1q','1l4j2b']:
                key=(ch,cat)
                try:
                    histos[proc][key]=fIn.Get(var+'_%s%s'%key).Clone(var+'_%s%s_'%key+proc)
                    histos[proc][key].SetDirectory(0)
                except:
                    pass
                if proc!='data' or cat!='1l4j2q': continue
                qcdcat=cat.replace('l','f')
                histos['qcd'][key]=fIn.Get('%s_%s%s'%(var,ch,qcdcat)).Clone('%s_%s%s_%s'%(var,ch,qcdcat,proc))
                histos['qcd'][key].SetDirectory(0)
        fIn.Close()

    #fix zeros
    for key in histos:
        for cat in histos[key]:
            histos[key][cat].Rebin()
            for xbin in xrange(1,histos[key][cat].GetNbinsX()+1):
                val=histos[key][cat].GetBinContent(xbin)
                histos[key][cat].SetBinContent(xbin,max(val,1e-3))

    return histos


"""
"""
def fitQCD(args,varName='mtw'):
    
    #fit the QCD yields
    histos=getHistos(args,varName)
    for key in histos['data']:

        var=ROOT.RooRealVar(varName,
                            'M_{T} [GeV]' if varName=='mtw' else 'MET [GeV]',
                            0,200)
        data=ROOT.RooDataHist('data','data',ROOT.RooArgList(var),histos['data'][key])

        ttbar=ROOT.RooDataHist('ttbar','ttbar',ROOT.RooArgList(var),histos['ttbar'][key])
        ttbarpdf=ROOT.RooHistPdf('ttbarpdf','ttbarpdf',ROOT.RooArgSet(var),ttbar,0)
        nttbar=ROOT.RooRealVar('nttbar','nttbar',0,1e6)

        wjetspdf=None
        try:
            wjets=ROOT.RooDataHist('wjets','wjets',ROOT.RooArgList(var),histos['wjets'][key])
            wjetspdf=ROOT.RooHistPdf('wjetspdf','wjetspdf',ROOT.RooArgSet(var),wjets,0)
            nwjets=ROOT.RooRealVar('nwjets','nwjets',0,1e6)
        except:
            pass

        qcdKey=(key[0],'1l4j2q')
        qcd=ROOT.RooDataHist('qcd','qcd',ROOT.RooArgList(var),histos['qcd'][qcdKey])
        qcdpdf=ROOT.RooHistPdf('qcdpdf','qcdpdf',ROOT.RooArgSet(var),qcd,0)
        nqcd=ROOT.RooRealVar('nqcd','nqcd',0,1e6)

        #sum up available contributions
        model=None
        if wjetspdf:
            model=ROOT.RooAddPdf('model','model',ROOT.RooArgList(ttbarpdf,wjetspdf,qcdpdf),ROOT.RooArgList(nttbar,nwjets,nqcd)) 
        else:
            model=ROOT.RooAddPdf('model','model',ROOT.RooArgList(ttbarpdf,qcdpdf),ROOT.RooArgList(nttbar,nqcd) )

        model.fitTo(data,ROOT.RooFit.Save(),ROOT.RooFit.Minos(1))

        #
        qcdEst,qcdEstUnc=nqcd.getVal(),nqcd.getError()
        if not key in QCDNORM:
            QCDNORM[key]=(qcdEst,qcdEstUnc)
        else:
            diff=abs(qcdEst-QCDNORM[key][0])
            QCDNORM[key]=(QCDNORM[key][0],
                          ROOT.TMath.Sqrt(diff**2+QCDNORM[key][1]**2))

        #show fit results
        c.Clear()
        frame=var.frame()
        data.plotOn(frame,ROOT.RooFit.Name('data'))
        color=ROOT.kGray
        model.plotOn(frame,
                     ROOT.RooFit.Components('tt*,w*'),
                     ROOT.RooFit.Name('ttw'),
                     ROOT.RooFit.ProjWData(data),
                     ROOT.RooFit.LineColor(color),
                     ROOT.RooFit.FillColor(color),
                     ROOT.RooFit.FillStyle(1001),
                     ROOT.RooFit.DrawOption('f'),
                     ROOT.RooFit.LineWidth(1),
                     ROOT.RooFit.MoveToBack()
                     )
        model.plotOn(frame,ROOT.RooFit.Name('qcd'),ROOT.RooFit.ProjWData(data),ROOT.RooFit.MoveToBack())
        frame.Draw()

        leg=ROOT.TLegend(0.6,0.7,0.9,0.6)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.AddEntry(data,'Data','ep')
        leg.AddEntry('qcd','QCD','l')
        leg.AddEntry('ttw','t#bar{t},W','f')
        leg.Draw()
    
        label = ROOT.TLatex()
        label.SetNDC()
        label.SetTextFont(42)
        label.SetTextSize(0.035)
        label.DrawLatex(0.6,0.9,'#bf{CMS} #it{preliminary}')
        label.DrawLatex(0.6,0.86,'174 nb^{-1} (#sqrt{s_{NN}}=8.16 TeV)')
        label.DrawLatex(0.6,0.82,'#it{%s%s}'%key)
        label.DrawLatex(0.6,0.78,'N_{QCD}=%3.0f#pm%3.0f'%(qcdEst,qcdEstUnc))

        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs(varName+'fit_%s%s.'%key+ext)
        #raw_input()


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.04)
c.SetBottomMargin(0.12)

fitQCD(sys.argv[1:],'met')
fitQCD(sys.argv[1:],'mtw')

#saturate the uncertainty at 50% at least but not more than 100%
for key in QCDNORM:
    val,valUnc=QCDNORM[key]
    valUnc=min(val,max(0.5*val,valUnc))
    QCDNORM[key]=(val,valUnc)

#dump results to a pickle file
with open('qcdnorm.pck','wb') as fOut:
    pickle.dump(QCDNORM,fOut)
