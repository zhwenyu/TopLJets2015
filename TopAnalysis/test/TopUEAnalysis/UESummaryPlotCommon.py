#!/usr/bin/env python

import ROOT
from UEAnalysisHandler import VARS,EVAXES

FILLS=[1001,3004,3002,1001]
COLORS=[ROOT.kAzure+4, ROOT.kMagenta, ROOT.kGreen+3,  ROOT.kMagenta+2, ROOT.kMagenta-9,ROOT.kRed+1,ROOT.kAzure+7, ROOT.kBlue-7,ROOT.kGray,ROOT.kGray]

MARKERS=[22,24,27,23,33,20,32,24]
OBSERVABLES=['chmult','sphericity','C','D','aplanarity','chavgpt','chavgpz','chflux','chfluxz']
OBSRANGES={'sphericity':(5e-2,5),
           'aplanarity':(5e-3,30),
           'C':(5e-2,5),
           'D':(5e-3,10),
           'chmult':(5e-4,1),
           'chavgpt':(2e-3,2),
           'chavgpz':(2e-3,2),
           'chflux':(5e-5,5e-2),
           'chfluxz':(5e-5,5e-2)}
RATIORANGES={'sphericity':(0.8,1.27),
           'aplanarity':(0.8,1.47),
           'C':(0.8,1.27),
           'D':(0.8,1.47),
           'chmult':(0.5,2.17),
           'chavgpt':(0.6,1.67),
           'chavgpz':(0.7,1.47),
           'chflux':(0.5,1.77),
           'chfluxz':(0.5,1.77)}

SLICES=['nj',None] #,'nj','ptttbar','ptll']
MAINMC=('POWHEG+PY8 CUETP8M2T4','t#bar{t}')
#MAINMC=('POWHEG+HW++ EE5C','t#bar{t} Herwig++')
COMPARISONSETS=[
    ('POWHEG+PY8 CUETP8M2T4', [ ('nominal',         ['t#bar{t}']), 
                                ('#deltaCUET8P2MT4',['t#bar{t} UEup',     't#bar{t} UEdn']),
                                ('FSR',             ['t#bar{t} fsr up',   't#bar{t} fsr dn']),
                                ('ISR',             ['t#bar{t} isr up',   't#bar{t} isr dn']),
                                ('hdamp',           ['t#bar{t} hdamp up', 't#bar{t} hdamp dn']),
                                ('CR',              ['t#bar{t} QCDbased', 't#bar{t} ERDon', 't#bar{t} gluon move']) ] 
     ),
    ('aMC@NLO+PY8 CUETP8M2T4', [ ('nominal', ['t#bar{t} aMC@NLO']) ]),
    ('POWHEG+HW++ EE5C'      , [ ('nominal', ['t#bar{t} Herwig++']) ]),
    ]

"""
"""
class SimpleUEPlot:
    def __init__(self,h,hname,trueAxis,addXunc):

        self.h=h.Clone(hname)
        self.h.SetDirectory(0)
        self.uncH={}
        self.trueAxis=trueAxis        
        self.addXunc=addXunc
        self.xmean=0
        self.x2mean=0
        self.generateGraph()
        self.garbage=[self.h,self.gr]

    def clear(self):
        for o in self.garbage:
            try:
                o.Delete()
            except:
                pass

    def subtractContribution(self,bkg):
        self.h.Add(bkg,-1)
        self.generateGraph()

    def generateGraph(self):
        self.gr=ROOT.TGraphAsymmErrors()
        self.gr.SetName('%s_gr'%self.h.GetName())

        norm=self.h.Integral()

        self.xmean=0
        self.x2mean=0
        for xbin in xrange(1,self.trueAxis.GetNbins()+1):
            cen,hwid=self.trueAxis.GetBinCenter(xbin),self.trueAxis.GetBinWidth(xbin)*0.5
            cts=self.h.GetBinContent(xbin)

            #include uncertainties (add envelopes in quadrature)
            errLo,errHi=self.h.GetBinError(xbin),self.h.GetBinError(xbin)
            for syst in self.uncH:
                minLo,maxHi=9999999999,-9999999999
                for k in xrange(0,len(self.uncH[syst])):
                    dSignal=expSysts[syst][k].GetBinContent(xbin)
                    if dSignal<minLo : minLo=dSignal
                    if dSignal>maxHi : maxHi=dSignal
                if minLo<0 and maxHi>0:
                    errLo += minLo**2
                    errHi += maxHi**2
                else:
                    env=max(abs(minLo),abs(maxHi))
                    errLo += env**2
                    errHi += env**2
            errLo=ROOT.TMath.Sqrt(errLo)
            errHi=ROOT.TMath.Sqrt(errHi)

            np=self.gr.GetN()
            self.gr.SetPoint(np, cen, cts/(norm*hwid) )
            self.gr.SetPointError(np,
                                  hwid if self.addXunc else 0,
                                  hwid if self.addXunc else 0,
                                  ROOT.TMath.Sqrt(errLo)/(norm*hwid),
                                  ROOT.TMath.Sqrt(errHi)/(norm*hwid))
    
"""
"""
class UESummaryPlotInfo:

    """
    """
    def __init__(self,obsAxis,sliceAxis):
        
        self.obsAxis=obsAxis
        self.sliceAxis=sliceAxis
        self.obs=obsAxis.GetName().split('_')[0]
        self.sliceVar=sliceAxis.GetName().split('_')[0] if sliceAxis else None

        self.data = []
        self.signal = []
        self.expSysts = []
        self.signalVars = []


    """
    """
    def addData(self,h,islice,ireg):
        self.data.append( (islice,ireg,
                           SimpleUEPlot(h,'data_%s_%s'%(str(islice),str(ireg)),self.obsAxis,False) ) )

    """
    """
    def addSignal(self,h,islice,ireg):
        self.signal.append( (islice,ireg,
                             SimpleUEPlot(h,'signal_%s_%s'%(str(islice),str(ireg)),self.obsAxis,True) ) )
        
    """
    subtracts the background
    """
    def subtractBackground(self,bkg,islice,ireg):
        for i,r,d in self.data:
            if i!=islice and r!=islice : continue
            d.subtractContribution(bkg)


    

    """
    """
    def show(self,outDir):

        sliceInfo={}
        for i,r,_ in self.data:
            if not i in sliceInfo: sliceInfo[i]=0
            sliceInfo[i]+=1

        for i in sliceInfo:
            npads=1+sliceInfo[i]

            cheight=200*(npads-1)+400
            c=ROOT.TCanvas('c','c',500,cheight)
            c.SetTopMargin(0.0)
            c.SetRightMargin(0.0)
            c.SetLeftMargin(0.0)
            c.SetBottomMargin(0.0)

            mheight=200.*(npads-1)/cheight
            print cheight,mheight,npads
            p1=ROOT.TPad('p1','p1',0.0,mheight,1.0,1.0) 
            p1.SetRightMargin(0.02)
            p1.SetLeftMargin(0.12)
            p1.SetTopMargin(0.06)
            p1.SetBottomMargin(0.01)
            p1.Draw()
            p1.cd()
            p1.SetLogy()
            frame=ROOT.TH1F('frame','frame',1,self.obsAxis.GetXmin(),self.obsAxis.GetXmax())
            frame.Draw()
            frame.GetYaxis().SetRangeUser(OBSRANGES[self.obs][0],OBSRANGES[self.obs][1])
            frame.GetYaxis().SetTitle('PDF')
            frame.GetYaxis().SetTitleOffset(1.0)
            frame.GetXaxis().SetLabelSize(0)
            frame.GetYaxis().SetTitleSize(0.05)
            frame.GetYaxis().SetLabelSize(0.05)

            mcleg=ROOT.TLegend(0.7,0.8,0.95,0.55)
            mcleg.SetFillStyle(0)
            mcleg.SetBorderSize(0)
            mcleg.SetTextFont(42)
            mcleg.SetTextSize(0.035)
            #draw the mominal prediction
            for j,r,p in self.signal:
                if i!=j : continue            

                p.gr.Draw('2')
                ci=ROOT.kAzure+7              
                regTitle=''
                if r==0: 
                    ci=ROOT.kMagenta
                    regTitle='toward'
                if r==1: 
                    ci=ROOT.kGreen+3
                    regTitle='transverse'
                if r==2: 
                    ci=ROOT.kAzure+4
                    regTitle='away'
                if r=='inc':
                    mcleg.AddEntry(p.gr,'#scale[0.8]{#splitline{Powheg+Pythia8}{CUETP8M2T4}}','f')
                else:
                    mcleg.AddEntry(p.gr,regTitle,'f')

                p.gr.SetMarkerStyle(1)
                p.gr.SetFillStyle(1001)
                p.gr.SetFillColor(ci)
                p.gr.SetLineColor(ci)
                p.gr.SetMarkerColor(ci)

            mcleg.Draw()

            #draw the data
            dleg=ROOT.TLegend(0.6,0.8,0.74,0.55)
            dleg.SetFillStyle(0)
            dleg.SetBorderSize(0)
            dleg.SetTextFont(42)
            dleg.SetTextSize(0.035)
            for j,r,p in self.data:
                if i!=j : continue
                p.gr.SetMarkerStyle(20)
                if r==0: 
                    p.gr.SetMarkerStyle(24)
                if r==1: 
                    p.gr.SetMarkerStyle(26)
                if r==2: 
                    p.gr.SetMarkerStyle(32)
                p.gr.Draw('ep')
                if r=='inc':
                    dleg.AddEntry(p.gr,'Data','ep')
                else:
                    dleg.AddEntry(p.gr,'','ep')
            dleg.Draw()

            #plot header
            tex=ROOT.TLatex()
            tex.SetTextFont(42)
            tex.SetTextSize(0.05)
            tex.SetNDC()
            tex.DrawLatex(0.7,0.85,'#bf{CMS} #it{preliminary}')
            #if opt.sliceVar:
            #    tex.DrawLatex(0.5,0.85,'%s #in %s'%(VARS[opt.sliceVar][0],idataGr.GetTitle()))
            tex.DrawLatex(0.74,0.955,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')

           
            rpads,rframes=[],[]
            for j in xrange(1,npads):
                c.cd()
                rheight=200./cheight
                rpads.append( ROOT.TPad('p2','p2',0.0,(j-1)*rheight,1.0,j*rheight) )
                rpads[-1].Draw()
                rpads[-1].SetRightMargin(0.02)
                rpads[-1].SetLeftMargin(0.12)
                rpads[-1].SetBottomMargin(0.2 if j==1 else 0.02)
                rpads[-1].SetTopMargin(0.01)
                rpads[-1].Draw()
                rpads[-1].cd()
                rframes.append( ROOT.TH1F('rframe%d'%j,'rframe%d'%j,1,self.obsAxis.GetXmin(),self.obsAxis.GetXmax()) )
                rframes[-1].GetYaxis().SetTitleSize(0.09)
                rframes[-1].GetYaxis().SetLabelSize(0.09)
                rframes[-1].GetYaxis().SetTitle('MC/Data')
                rframes[-1].GetYaxis().SetTitleOffset(0.7)
                rframes[-1].GetXaxis().SetTitleSize(0.09 if j==1 else 0.0)
                rframes[-1].GetXaxis().SetLabelSize(0.09 if j==1 else 0.0)
                rframes[-1].GetXaxis().SetTitle(VARS[self.obs][0])
                rframes[-1].GetYaxis().SetRangeUser(RATIORANGES[self.obs][0],RATIORANGES[self.obs][1])
                rframes[-1].GetYaxis().SetNdivisions(5+100*5)
                rframes[-1].Draw()


            c.cd()
            c.Modified()
            c.Update()
            raw_input()


#
#
#
#    for islice in xrange(1,nslices+1):
#
#        idataGr=ROOT.TGraphErrors();
#        idataGr.SetMarkerStyle(20+(islice-1))
#        idataGr.SetMarkerColor(1)
#        idataGr.SetLineColor(1)
#        idataGr.SetName('data_%d'%islice)
#        title='inc'
#        if sliceAxis:
#            title='[%d,%d]'%(int(sliceAxis.GetBinLowEdge(islice)),int(sliceAxis.GetBinUpEdge(islice)))
#        idataGr.SetTitle(title)
#
#        isignalGr=ROOT.TGraphAsymmErrors();
#        isignalGr.SetFillStyle(1001)
#        ci=ROOT.kAzure+7
#        isignalGr.SetFillColor(ci)
#        isignalGr.SetMarkerColor(ci)
#        isignalGr.SetLineColor(ci)
#        isignalGr.SetMarkerStyle(1)
#
#        isignalRatioGr=isignalGr.Clone('nominalratio')
#
#
#        ratiosGr=ROOT.TMultiGraph()
#        #ratiosLeg=ROOT.TLegend(0.7,0.2,0.95,0.6)
#        ratiosLeg=ROOT.TLegend(0.1,0.85,0.98,0.9)
#        ratiosLeg.SetNColumns( len(signalVars)+1 )
#
#        #build final distribution
#        nobsBins=obsAxis.GetNbins()
#        for xbin in xrange(1,nobsBins+1):
#            cen,hwid=obsAxis.GetBinCenter(xbin),obsAxis.GetBinWidth(xbin)*0.5
#            
#            rawBin=xbin
#            if sliceAxis: rawBin+=nobsBins*(islice-1)
#            
#            #experimental systematics (add envelopes in quadrature)
#            errLo,errHi=0,0
#            for syst in expSysts:
#                minLo,maxHi=9999999999,-9999999999
#                for k in xrange(0,len(expSysts[syst])):
#                    dSignal=expSysts[syst][k].GetBinContent(rawBin)
#                    if dSignal<minLo : minLo=dSignal
#                    if dSignal>maxHi : maxHi=dSignal
#                if minLo<0 and maxHi>0:
#                    errLo += minLo**2
#                    errHi += maxHi**2
#                else:
#                    env=max(abs(minLo),abs(maxHi))
#                    errLo += env**2
#                    errHi += env**2
#            errLo=ROOT.TMath.Sqrt(errLo)
#            errHi=ROOT.TMath.Sqrt(errHi)
#
#            np=idataGr.GetN()
#
#            dataCts,dataUnc=data.GetBinContent(rawBin),data.GetBinError(rawBin)
#            idataGr.SetPoint       (np,   cen,        dataCts/(2*hwid))
#            idataGr.SetPointError  (np,   0,          dataUnc/(2*hwid))
#
#            signalCts=signal.GetBinContent(rawBin)
#            isignalGr.SetPoint     (np,   cen,        signalCts/(2*hwid))
#            isignalGr.SetPointError(np,   hwid, hwid, errLo/(2*hwid), errHi/(2*hwid))
#
#            if dataCts>0 and signalCts>0:
#                iratioVal=signalCts/dataCts
#                iratioUncLo=iratioVal*ROOT.TMath.Sqrt(pow(errLo/signalCts,2)+pow(dataUnc/dataCts,2))
#                iratioUncHi=iratioVal*ROOT.TMath.Sqrt(pow(errHi/signalCts,2)+pow(dataUnc/dataCts,2))
#                isignalRatioGr.SetPoint     (np,   cen,        iratioVal)
#                isignalRatioGr.SetPointError(np,   hwid, hwid, iratioUncLo, iratioUncHi)
#
#        ratiosGr.Add(isignalRatioGr,'2')
#        ratiosLeg.AddEntry(isignalRatioGr,'PW+PY8 CUETP8M2T4','f')
#
#        #dataAvg=averageDistribution(idataGr,obsAxis)
#        #mcAvg={'nominal':averageDistribution(isignalGr,obsAxis)}
#
#        #variations to be compared
#        ivar=0
#        for var in signalVars:
#            ivar+=1
#            isignalVarGr=ROOT.TGraphAsymmErrors();
#            isignalVarGr.SetFillStyle(1001)
#            ci=ROOT.TColor.GetColor(COLORS[ivar-1])
#            isignalVarGr.SetFillColor(ci)
#            isignalVarGr.SetMarkerColor(ci)
#            isignalVarGr.SetLineColor(ci)
#            isignalVarGr.SetMarkerStyle(MARKERS[ivar-1])
#            isignalVarGr.SetMarkerSize(0.8)
#
#            for xbin in xrange(1,nobsBins+1):
#                cen,hwid=obsAxis.GetBinLowEdge(xbin),obsAxis.GetBinWidth(xbin)*0.5
#                cen += 2*hwid*ivar/(len(signalVars)+2)
#
#                rawBin=xbin
#                if sliceAxis: rawBin+=nobsBins*(islice-1)
#
#                minR2data,maxR2data=200,-200
#                for varH in signalVars[var]:
#                    r2data=varH.GetBinContent(rawBin)
#                    minR2data=min(r2data,minR2data)
#                    maxR2data=max(r2data,maxR2data)
#                
#                y=0.5*(minR2data+maxR2data)
#                errHi=max(maxR2data-y,1e-4) #just for display purposes
#                errLo=max(y-minR2data,1e-4)
#                            
#                np=isignalVarGr.GetN()
#                isignalVarGr.SetPoint(np,cen,y)
#                isignalVarGr.SetPointError(np,0,0,errLo,errHi)
#            
#            ratiosGr.Add(isignalVarGr,'pZ')
#            ratiosLeg.AddEntry(isignalVarGr,var,'ep')
#
#
#        

#        ratiosGr.Draw('p')
#        l=ROOT.TLine(frameratio.GetXaxis().GetXmin(),1,frameratio.GetXaxis().GetXmax(),1)
#        l.SetLineColor(ROOT.kBlue)
#        l.Draw()
#        ratiosLeg.SetBorderSize(0)
#        ratiosLeg.SetFillStyle(0)
#        ratiosLeg.SetTextFont(42)
#        ratiosLeg.SetTextSize(0.06)
#        ratiosLeg.Draw()
#        
#        c.cd()
#        c.Modified()
#        c.Update()
#        outname=obs
#        if nslices>1: outname += '%s_%d'%(sliceVar,islice)
#        for ext in ['png','pdf']:
#            c.SaveAs('~/www/TopUE_ReReco2016/%s_smeared.%s'%(outname,ext))
#

