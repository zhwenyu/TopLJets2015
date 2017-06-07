import ROOT
import os
import sys
import optparse
import pickle
from collections import OrderedDict
from UEAnalysisHandler import VARS

COLORS=[ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kMagenta-9,ROOT.kRed+1,ROOT.kAzure+7, ROOT.kBlue-7,ROOT.kGray,ROOT.kGray]
MARKERS=[22,24,20,27,23,33,20,32,24]
OBSERVABLES=['chmult','sphericity','C','D','aplanarity','chavgpt','chavgpz','chflux','chfluxz']
OBSRANGES={'sphericity':(5e-3,5),
           'aplanarity':(5e-3,30),
           'C':(5e-3,5),
           'D':(5e-3,10),
           'chmult':(5e-4,0.1),
           'chavgpt':(2e-3,2),
           'chavgpz':(2e-3,2),
           'chflux':(5e-5,5e-2),
           'chfluxz':(5e-5,5e-2)}
RATIORANGES={'sphericity':(0.8,1.27),
           'aplanarity':(0.8,1.27),
           'C':(0.8,1.27),
           'D':(0.8,1.27),
           'chmult':(0.5,2.17),
           'chavgpt':(0.5,1.47),
           'chavgpz':(0.5,1.47),
           'chflux':(0.5,1.97),
           'chfluxz':(0.5,1.87)}

SLICES=[None,'nj','ptttbar'] #,'ptll']

"""
"""
def averageDistribution(gr,axis):
    x,y=ROOT.Double(0),ROOT.Double(0)
    avgX,avgXhi,avgXlo=0
    totaly=0
    for ip in xrange(0,gr.GetN()):
           gr.GetPoint(ip,x,y)       
           x=axis.GetBinCenter(ip+1)
           wid=axis.GetBinWidth(ip+1)
           eyhi,eylo=gr.GetEYhigh(),gr.GetEYlow()
           totaly += y*wid
           avgX   += x*y*wid
           avgXhi += x*(y+eyhi)*wid
           avgXlo += x*(y-eylo)*wid
    return avgX/totaly,avgXhi/totaly,avgXlo/totaly
    

"""
"""
def buildPlot(data,signal,expSysts,signalVars,obsAxis,sliceAxis,opt):

    obs=obsAxis.GetName().split('_')[0]
    sliceVar=sliceAxis.GetName().split('_')[0] if sliceAxis else ''
    nslices=sliceAxis.GetNbins() if sliceAxis else 1
    frame=ROOT.TH1F('frame','frame',1,obsAxis.GetXmin(),obsAxis.GetXmax())
    frameratio=ROOT.TH1F('frameratio','frameratio',1,obsAxis.GetXmin(),obsAxis.GetXmax())

    for islice in xrange(1,nslices+1):

        idataGr=ROOT.TGraphErrors();
        idataGr.SetMarkerStyle(20+(islice-1))
        idataGr.SetMarkerColor(1)
        idataGr.SetLineColor(1)
        idataGr.SetName('data_%d'%islice)
        title='inc'
        if sliceAxis:
            title='[%d,%d]'%(int(sliceAxis.GetBinLowEdge(islice)),int(sliceAxis.GetBinUpEdge(islice)))
        idataGr.SetTitle(title)

        isignalGr=ROOT.TGraphAsymmErrors();
        isignalGr.SetFillStyle(1001)
        ci=ROOT.kAzure+7
        isignalGr.SetFillColor(ci)
        isignalGr.SetMarkerColor(ci)
        isignalGr.SetLineColor(ci)
        isignalGr.SetMarkerStyle(1)

        isignalRatioGr=isignalGr.Clone('nominalratio')
        
        ratiosGr=ROOT.TMultiGraph()
        #ratiosLeg=ROOT.TLegend(0.7,0.2,0.95,0.6)
        ratiosLeg=ROOT.TLegend(0.1,0.85,0.98,0.9)
        ratiosLeg.SetNColumns( len(signalVars)+1 )

        #build final distribution
        nobsBins=obsAxis.GetNbins()
        for xbin in xrange(1,nobsBins+1):
            cen,hwid=obsAxis.GetBinCenter(xbin),obsAxis.GetBinWidth(xbin)*0.5
            
            rawBin=xbin
            if sliceAxis: rawBin+=nobsBins*(islice-1)
            
            #experimental systematics (add envelopes in quadrature)
            errLo,errHi=0,0
            for syst in expSysts:
                minLo,maxHi=9999999999,-9999999999
                for k in xrange(0,len(expSysts[syst])):
                    dSignal=expSysts[syst][k].GetBinContent(rawBin)
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

            np=idataGr.GetN()

            dataCts,dataUnc=data.GetBinContent(rawBin),data.GetBinError(rawBin)
            idataGr.SetPoint       (np,   cen,        dataCts/(2*hwid))
            idataGr.SetPointError  (np,   0,          dataUnc/(2*hwid))

            signalCts=signal.GetBinContent(rawBin)
            isignalGr.SetPoint     (np,   cen,        signalCts/(2*hwid))
            isignalGr.SetPointError(np,   hwid, hwid, errLo/(2*hwid), errHi/(2*hwid))

            if dataCts>0 and signalCts>0:
                iratioVal=signalCts/dataCts
                iratioUncLo=iratioVal*ROOT.TMath.Sqrt(pow(errLo/signalCts,2)+pow(dataUnc/dataCts,2))
                iratioUncHi=iratioVal*ROOT.TMath.Sqrt(pow(errHi/signalCts,2)+pow(dataUnc/dataCts,2))
                isignalRatioGr.SetPoint     (np,   cen,        iratioVal)
                isignalRatioGr.SetPointError(np,   hwid, hwid, iratioUncLo, iratioUncHi)


        ratiosGr.Add(isignalRatioGr,'2')
        ratiosLeg.AddEntry(isignalRatioGr,'PW+PY8 CUETP8M2T4','f')

        dataAvg=averageDistribution(idataGr,obsAxis)
        mcAvg={'nominal':averageDistribution(isignalGr,obsAxis)}

        #variations to be compared
        ivar=0
        for var in signalVars:
            ivar+=1
            isignalVarGr=ROOT.TGraphAsymmErrors();
            isignalVarGr.SetFillStyle(1001)
            ci=ROOT.TColor.GetColor(COLORS[ivar-1])
            isignalVarGr.SetFillColor(ci)
            isignalVarGr.SetMarkerColor(ci)
            isignalVarGr.SetLineColor(ci)
            isignalVarGr.SetMarkerStyle(MARKERS[ivar-1])
            isignalVarGr.SetMarkerSize(0.8)

            for xbin in xrange(1,nobsBins+1):
                cen,hwid=obsAxis.GetBinLowEdge(xbin),obsAxis.GetBinWidth(xbin)*0.5
                cen += 2*hwid*ivar/(len(signalVars)+2)

                rawBin=xbin
                if sliceAxis: rawBin+=nobsBins*(islice-1)

                minR2data,maxR2data=200,-200
                for varH in signalVars[var]:
                    r2data=varH.GetBinContent(rawBin)
                    minR2data=min(r2data,minR2data)
                    maxR2data=max(r2data,maxR2data)
                
                y=0.5*(minR2data+maxR2data)
                errHi=max(maxR2data-y,1e-4) #just for display purposes
                errLo=max(y-minR2data,1e-4)
                            
                np=isignalVarGr.GetN()
                isignalVarGr.SetPoint(np,cen,y)
                isignalVarGr.SetPointError(np,0,0,errLo,errHi)
            
            ratiosGr.Add(isignalVarGr,'pZ')
            ratiosLeg.AddEntry(isignalVarGr,var,'ep')

        c=ROOT.TCanvas('c','c',1000,600)
        c.SetTopMargin(0.0)
        c.SetRightMargin(0.0)
        c.SetLeftMargin(0.0)
        c.SetBottomMargin(0.0)
        p1=ROOT.TPad('p1','p1',0.0,0.4,1.0,1.0) 
        p1.SetRightMargin(0.02)
        p1.SetLeftMargin(0.07)
        p1.SetTopMargin(0.05)
        p1.SetBottomMargin(0.01)
        p1.Draw()
        p1.cd()
        p1.SetLogy()
        frame.Draw()
        frame.GetYaxis().SetRangeUser(OBSRANGES[obs][0],OBSRANGES[obs][1])
        frame.GetYaxis().SetTitle('PDF')
        frame.GetYaxis().SetTitleOffset(0.6)
        frame.GetXaxis().SetLabelSize(0)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.05)
        isignalGr.Draw('2')
        idataGr.Draw('ep')
        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.07)
        tex.SetNDC()
        tex.DrawLatex(0.1,0.85,'#bf{CMS} #it{preliminary}')
        if nslices>1:
            tex.DrawLatex(0.5,0.85,'%s #in %s'%(VARS[sliceVar][0],idataGr.GetTitle()))
        tex.DrawLatex(0.85,0.955,'#scale[0.6]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')
        leg=ROOT.TLegend(0.7,0.9,0.9,0.65)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.05)
        leg.AddEntry(idataGr,'Data','ep')
        leg.AddEntry(isignalGr,'#splitline{Powheg+Pythia8}{CUETP8M2T4}','f')
        leg.Draw()
        
        c.cd()
        p2 = ROOT.TPad('p2','p2',0.0,0.0,1.0,0.4)
        p2.Draw()
        p2.SetBottomMargin(0.18)
        p2.SetRightMargin(0.02)
        p2.SetLeftMargin(0.07)
        p2.SetTopMargin(0.01)
        p2.Draw()
        p2.cd()
        frameratio.Draw()
        frameratio.GetYaxis().SetTitleSize(0.07)
        frameratio.GetYaxis().SetLabelSize(0.07)
        frameratio.GetXaxis().SetTitleSize(0.07)
        frameratio.GetXaxis().SetLabelSize(0.07)
        frameratio.GetYaxis().SetTitle('MC/Data')
        frameratio.GetYaxis().SetTitleOffset(0.5)
        frameratio.GetXaxis().SetTitle(VARS[obs][0])
        frameratio.GetYaxis().SetRangeUser(RATIORANGES[obs][0],RATIORANGES[obs][1])
        frameratio.GetYaxis().SetNdivisions(5+100*5)
        ratiosGr.Draw('p')
        l=ROOT.TLine(frameratio.GetXaxis().GetXmin(),1,frameratio.GetXaxis().GetXmax(),1)
        l.SetLineColor(ROOT.kBlue)
        l.Draw()
        ratiosLeg.SetBorderSize(0)
        ratiosLeg.SetFillStyle(0)
        ratiosLeg.SetTextFont(42)
        ratiosLeg.SetTextSize(0.06)
        ratiosLeg.Draw()
        
        c.cd()
        c.Modified()
        c.Update()
        outname=obs
        if nslices>1: outname += '%s_%d'%(sliceVar,islice)
        c.SaveAs('~/www/TopUE_ReReco2016/%s.png'%outname)



def normalizePerSlice(h,obsAxis,sliceAxis):

    nobsBins=obsAxis.GetNbins()
    nslices=sliceAxis.GetNbins() if sliceAxis else 1
    norm=[0]*nslices
    for islice in xrange(1,nslices+1):
        for xbin in xrange(1,nobsBins+1):
            rawBin=xbin
            if sliceAxis: rawBin+=nobsBins*(islice-1)
            norm[islice-1]+=h.GetBinContent(rawBin)
    for islice in xrange(1,nslices+1):
        for xbin in xrange(1,nobsBins+1):
            rawBin=xbin
            if sliceAxis: rawBin+=nobsBins*(islice-1)
            h.SetBinContent( rawBin, h.GetBinContent(rawBin)/norm[islice-1] )
            h.SetBinError( rawBin, h.GetBinError(rawBin)/norm[islice-1] )

"""
"""
def readPlotsFrom(args,opt):

    outdir=args[0].replace('.root','')

    #analysis axes
    analysisaxis=None
    with open(opt.analysisAxis,'r') as cachefile:
        analysisaxis = pickle.load(cachefile)

    #variants to compare to
    varTypes=OrderedDict()
    for varDesc in opt.vars.split(';'):
        try:
            varTitle,varNamesEnum=varDesc.split(':')
            varList=[x for x in varNamesEnum.split(',')]
            varTypes[varTitle]=varList
        except:
            pass

    fIn=ROOT.TFile.Open(args[0])
    fSyst=ROOT.TFile.Open(args[1]) if len(args)>1 else None
    for obs in OBSERVABLES:

        obsAxis=analysisaxis[(obs,True)]
        
        for s in SLICES:

            sliceAxis=None if s is None else analysisaxis[(s,True)]

            #read the nominal expectations
            nomKey='%s_%s_inc_None_True'%(obs,s)
            t=fIn.Get(nomKey)
            print t,nomKey
            data,signal,bkg=None,None,None
            for pkey in t.GetListOfKeys():
                h=t.Get(pkey.GetName())
                if not h.InheritsFrom('TH1') : continue
                if 'Data' in h.GetTitle():
                    data=h.Clone('data')
                elif h.GetTitle() in opt.signal:
                    signal=h.Clone('%s_nominal'%h.GetName())
                else:
                    if bkg is None: bkg=h.Clone('bkg')
                    else : bkg.Add(h)

            #normalize the signal
            normalizePerSlice(signal,obsAxis,sliceAxis)

            #subtract the background from the data 
            data.Add(bkg,-1)            
            normalizePerSlice(data,obsAxis,sliceAxis)


            #project experimental systematics and signal variations
            signalVars={}
            expSysts={}
            expSystsKey='%s_%s_inc_syst_True'%(obs,s)
            expSystsH=fIn.Get('%s/%s_%s'%(expSystsKey,expSystsKey,opt.signal))
            for ybin in xrange(2,expSystsH.GetNbinsY()):
                varName=expSystsH.GetYaxis().GetBinLabel(ybin)
                systKey=varName
                if systKey[-2:] in ['up','dn']  : systKey=systKey[:-2]
                h=expSystsH.ProjectionX('px',ybin,ybin)
                normalizePerSlice(h,obsAxis,sliceAxis)
                if systKey in ['mur','muf','q'] :                    
                    h.Divide(data)
                    systKey='ME scale'
                    if not systKey in signalVars: signalVars[systKey]=[]
                    signalVars[systKey].append( h.Clone(varName) )
                elif systKey in ['p_{T}(t)']:
                    h.Divide(data)
                    systKey='toppt'
                    if not systKey in signalVars: signalVars[systKey]=[]
                    signalVars[systKey].append( h.Clone(varName) )                    
                elif systKey in ['tkeff','tkeffbcdef','tkeffgh','tkeffeta']:
                    if systKey in ['tkeffbcdef','tkeffgh'] : continue
                    h.Add(signal,-1)
                    systKey='Trk. eff.'
                    if not systKey in expSysts: expSysts[systKey]=[]
                    expSysts[systKey].append( h.Clone(varName) )
                else:
                    h.Add(signal,-1)
                    if not systKey in expSysts: expSysts[systKey]=[]
                    expSysts[systKey].append(  h.Clone(varName) )

            #variations to compare to
            for varTitle in varTypes:
                signalVars[varTitle]=[]
                for varName in varTypes[varTitle]:
                    hvar=fSyst.Get(nomKey).Get(signal.GetName().replace('_nominal',' '+varName))
                    normalizePerSlice(hvar,obsAxis,sliceAxis)
                    hvar.Divide(data)
                    signalVars[varTitle].append(hvar)


            buildPlot(data,signal,expSysts,signalVars,obsAxis,sliceAxis,opt)




"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--signal',
                      dest='signal',
                      help='signal [%default]',
                      type='string',
                      default='t#bar{t}')
    parser.add_option('--vars',
                      dest='vars',
                      help='variations to outsource from other files [%default]',
                      default='#deltaCUET8P2MT4:UEup,UEdn;FSR:fsr up,fsr dn;ISR:isr up,isr dn;CR:QCDbased,ERDon;hdamp:hdamp up,hdamp dn;', #HW++EE5C:Herwig++;aMC@NLO:aMC@NLO',
                      type='string')
    parser.add_option('--cfg',
                      dest='analysisAxis',
                      help='cfg with axis definitions [%default]',
                      type='string',
                      default='%s/src/TopLJets2015/TopAnalysis/UEanalysis/analysisaxiscfg.pck'%os.environ['CMSSW_BASE'])
    (opt, args) = parser.parse_args()


    readPlotsFrom(args,opt)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
