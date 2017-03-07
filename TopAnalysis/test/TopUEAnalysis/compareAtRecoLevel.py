import ROOT
import os
import sys
import optparse
import numpy as np
from collections import OrderedDict

COLORS=['#d73027','#00374A','#fee090']
FILLS =[3001,     3001,     3001]

"""
returns a plot with the relative uncertainty
"""
def getRelUncPlot(h,name,title='rel unc',color=ROOT.TColor.GetColor('#99d8c9'),fill=1001):
    relUncH=h.Clone(name)
    relUncH.Reset('ICE')
    relUncH.SetDirectory(0)
    relUncH.SetTitle(title)
    relUncH.SetLineColor(color)
    relUncH.SetFillColor(color)
    relUncH.SetMarkerColor(color)
    relUncH.SetMarkerStyle(0)
    relUncH.SetFillStyle(fill)
    for xbin in xrange(1,h.GetNbinsX()+1):
        relUncH.SetBinContent(xbin,1)
        val=h.GetBinContent(xbin)
        unc=h.GetBinError(xbin)
        if val==0 : continue
        relUncH.SetBinError(xbin,unc/val)
    return relUncH

"""
"""
def getRatioPlotsFrom(pset1,pset2):
    allPlots={}
    for key in pset1:
        if not key in pset2: continue
        allPlots[key]=VarPlot(key,pset1[key].varH)
        nomRatio=pset1[key].nomH.Clone()
        nomRatio.Divide(pset2[key].nomH)
        nomRatio.GetYaxis().SetTitle('Ratio')
        allPlots[key].addNominal(nomRatio)

        allRatios=[]
        for varTitle in pset1[key].varH:
            if not varTitle in pset2[key].varH : continue
            allPlots[key].addCorrelation(varTitle)
            for var in pset1[key].varH[varTitle]:
                if not var in  pset2[key].varH[varTitle] : continue
                allRatios.append( pset1[key].varH[varTitle][var].Clone() )
                allRatios[-1].Divide( pset2[key].varH[varTitle][var] )
                allRatios[-1].GetYaxis().SetTitle('Ratio')
                allPlots[key].addVariation( allRatios[-1],var,varTitle )

                for xbin in xrange(1,pset1[key].varH[varTitle][var].GetNbinsX()):
                    xref=pset1[key].nomH.GetBinContent(xbin)
                    yref=pset2[key].nomH.GetBinContent(xbin)
                    if xref==0 or yref==0 : continue
                    xcen=pset1[key].varH[varTitle][var].GetBinContent(xbin)/xref
                    xunc=0.5*pset1[key].varH[varTitle][var].GetBinError(xbin)/xref
                    ycen=pset2[key].varH[varTitle][var].GetBinContent(xbin)/yref
                    yunc=0.5*pset2[key].varH[varTitle][var].GetBinError(xbin)/yref
                    allPlots[key].addPointToCorrelation(varTitle,xcen,xunc,ycen,yunc)

    return allPlots

"""
"""
class VarPlot:
    def __init__(self,var,varKeys):
        self.var=var
        self.nomH=None
        self.statUncRatioH=None
        self.varH=OrderedDict()
        self.corrH=OrderedDict()
        for key in varKeys : self.varH[key]=OrderedDict()
    def addNominal(self,h):
        self.nomH=h.Clone('nominal_%s'%self.var)
        self.nomH.SetDirectory(0)
        self.statUncRatioH=getRelUncPlot(self.nomH,
                                         name='relstatunc_%s'%self.var,
                                         title='Stat.',
                                         color=ROOT.TColor.GetColor('#99d8c9'),
                                         fill=1001)

    def getRelativeVariations(self):

        relVars=[]
        for varTitle in self.varH:
            irelVars=[]
            for var in self.varH[varTitle]:
                hname=self.varH[varTitle][var].GetName()
                irelVars.append( self.varH[varTitle][var].Clone('rel_%s'%hname) )
                #irelVars[-1].Divide(self.nomH)
                irelVars[-1].SetDirectory(0)

            #start a new graph
            iv=len(relVars)
            relVars.append( ROOT.TGraphErrors() )
            relVars[-1].SetTitle(varTitle)
            relVars[-1].SetFillColor(ROOT.TColor.GetColor(COLORS[iv%2]))
            relVars[-1].SetLineColor(ROOT.TColor.GetColor(COLORS[iv%2]))
            relVars[-1].SetMarkerColor(ROOT.TColor.GetColor(COLORS[iv%2]))
            relVars[-1].SetFillStyle(FILLS[iv%2])

            #fill the graph
            if len(irelVars)==2:

                hmean=irelVars[0].Clone('hmean')
                hmean.Add(irelVars[1])
                hmean.Scale(0.5)
                hdiff=irelVars[0].Clone('hdiff')
                hdiff.Add(irelVars[1],-1)
                hdiff.Scale(0.5)

                for xbin in xrange(1,hmean.GetNbinsX()+1):
                    x,dx=hmean.GetXaxis().GetBinCenter(xbin),0.5*hmean.GetXaxis().GetBinWidth(xbin)
                    y,dy=hmean.GetBinContent(xbin),ROOT.TMath.Abs(hdiff.GetBinContent(xbin))
                    nomy=self.nomH.GetBinContent(xbin)
                    relVars[-1].SetPoint(xbin-1,x,y/nomy)
                    relVars[-1].SetPointError(xbin-1,dx,dy/nomy)

                hmean.Delete()
                hdiff.Delete()
            elif len(irelVars)==1:
                relVars[-1].SetLineWidth(2)
                for xbin in xrange(1,irelVars[0].GetNbinsX()+1):
                    x,dx=irelVars[0].GetXaxis().GetBinCenter(xbin),0.5*irelVars[0].GetXaxis().GetBinWidth(xbin)
                    y=irelVars[0].GetBinContent(xbin)
                    nomy=self.nomH.GetBinContent(xbin)
                    relVars[-1].SetPoint(xbin-1,x,y/nomy)
                    relVars[-1].SetPointError(xbin-1,dx,0)
            elif 'NNPDF' in varTitle:

                hmean=self.nomH.Clone('hmean')

                #the last two are treated as the alphaS variations
                #hmean=irelVars[-1].Clone('hmean')
                #hmean.Add(irelVars[-2])
                #hmean.Scale(0.5)
                hdiff=irelVars[-1].Clone('hdiff')
                hdiff.Add(irelVars[-2],-1)
                hdiff.Scale(0.5)

                #all the others compute standard deviation and add in quadrature to the alphaS variation
                npdfvars=len(relVars)
                for xbin in xrange(1,irelVars[0].GetNbinsX()+1):


                    pdfVarVals=[]
                    for i in xrange(0,npdfvars-2) : pdfVarVals.append( irelVars[i].GetBinContent(xbin))
                    rms2=np.var(pdfVarVals)

                    x,dx=hmean.GetXaxis().GetBinCenter(xbin),0.5*hmean.GetXaxis().GetBinWidth(xbin)
                    y,dy=hmean.GetBinContent(xbin),ROOT.TMath.Sqrt(hdiff.GetBinContent(xbin)**2+rms2)
                    nomy=self.nomH.GetBinContent(xbin)

                    relVars[-1].SetPoint(xbin-1,x,y/nomy)
                    relVars[-1].SetPointError(xbin-1,dx,dy/nomy)

                hmean.Delete()
                hdiff.Delete()

        return relVars

    def addVariation(self,h,var,varTitle):
        if not varTitle in self.varH: self.varH[varTitle]={}
        self.varH[varTitle][var]=h.Clone('variation_%s_%s'%(self.var,var))
        self.varH[varTitle][var].SetDirectory(0)
    def Print(self):
        print '%s plot'%self.var
        if self.nomH : print ' contains nominal distribution'
        print ' contains %d variations'%len(self.varH)
    def show(self,outdir):

        if not self.nomH : return

        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)

        garbageList=[]

        #prepare output
        os.system('mkdir -p %s'%outdir)

        nVarSets=int((len(self.varH)-1)/2)+1

        #start canvas
        canvasY=200+300*nVarSets
        c = ROOT.TCanvas('c','c',500,canvasY)
        c.SetBottomMargin(0.0)
        c.SetLeftMargin(0.0)
        c.SetTopMargin(0)
        c.SetRightMargin(0.0)
        c.cd()
        garbageList.append(c)

        deltaY=0.08
        paddy=((canvasY*(1-deltaY))/float(nVarSets+1))/canvasY

        subPads,allLegs=[],[]


        #main plot
        p1=ROOT.TPad('p1','p1',0.0,1.0-paddy,1.0,1.0)
        p1.Draw()
        p1.SetRightMargin(0.05)
        p1.SetLeftMargin(0.12)
        p1.SetTopMargin(0.03)
        p1.SetBottomMargin(0.03)
        p1.SetGridx(False)
        p1.SetGridy(False)
        p1.cd()
        garbageList.append(p1)

        #frame
        frame = self.nomH.Clone('frame')
        frame.Reset('ICE')
        frame.SetDirectory(0)
        frame.GetYaxis().SetTitle('1/N dN/dx')
        frame.GetYaxis().SetTitleSize(0.1)
        frame.GetYaxis().SetLabelSize(0.09)
        frame.GetYaxis().SetNoExponent()
        frame.GetYaxis().SetTitleOffset(0.6)
        frame.GetXaxis().SetTitleSize(0.0)
        frame.GetXaxis().SetLabelSize(0.0)
        frame.Draw()
        garbageList.append(frame)

        #label
        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(42)
        txt.SetTextSize(0.1)
        txt.SetTextAlign(12)
        txt.SetTextSize(0.11)
        txt.DrawLatex(0.6,0.9,'#bf{CMS} #it{preliminary}')
        txt.DrawLatex(0.6,0.8,'#it{reconstruction level}')
        txt.DrawLatex(0.6,0.7,'37 fb^{-1} (13 TeV)')
        allLegs.append( ROOT.TLegend(0.6,0.65,0.9,0.4) )
        allLegs[-1].SetBorderSize(0)
        allLegs[-1].SetFillStyle(0)
        allLegs[-1].SetTextFont(42)
        allLegs[-1].SetTextSize(0.11)
        allLegs[-1].AddEntry(self.nomH,self.nomH.GetTitle(),'p')

        #nominal histogram
        firstVar=self.varH.keys()[0]
        for var in self.varH[firstVar]:
            hname=self.varH[firstVar][var].Draw('histsame')
            allLegs[-1].AddEntry(self.varH[firstVar][var],self.varH[firstVar][var].GetTitle(),'l')

        self.nomH.SetLineColor(1)
        self.nomH.SetLineWidth(2)
        self.nomH.Draw('e1same')
        frame.GetYaxis().SetRangeUser(self.nomH.GetMinimum()*0.9,
                                      self.nomH.GetMaximum()*1.2)

        allLegs[-1].Draw()

        #ratios
        ratioframe=frame.Clone('ratioframe')
        ratioframe.GetYaxis().SetTitle('Data/MC')
        ratioframe.GetXaxis().SetTitle(self.var)
        ratioframe.GetYaxis().SetRangeUser(0.65,1.55)
        ratioframe.GetYaxis().SetNdivisions(5)
        ratioframe.GetYaxis().SetLabelSize(0.09)
        ratioframe.GetYaxis().SetTitleSize(0.1)
        ratioframe.GetYaxis().SetTitleOffset(0.55)
        ratioframe.GetXaxis().SetLabelSize(0.09)
        ratioframe.GetXaxis().SetTitleSize(0.1)
        ratioframe.GetXaxis().SetTitleOffset(1.1)
        ratioframe.GetYaxis().SetNoExponent()
        ratioframe.SetFillStyle(3001)
        ratioframe.SetFillColor(ROOT.kGray+2)
        garbageList.append(ratioframe)

        relVars=self.getRelativeVariations()

        for ivarSet in xrange(0,nVarSets):
            c.cd()

            if ivarSet==nVarSets-1:
                subPads.append( ROOT.TPad('p2%d'%ivarSet,'p2%d'%ivarSet,0.0,0.0,1.0,1.0-paddy*(ivarSet+1)) )
            else:
                subPads.append( ROOT.TPad('p2%d'%ivarSet,'p2%d'%ivarSet,0.0,1.0-paddy*(ivarSet+2),1.0,1.0-paddy*(ivarSet+1)) )

            subPads[-1].Draw()
            subPads[-1].SetRightMargin(0.05)
            subPads[-1].SetLeftMargin(0.12)
            subPads[-1].SetTopMargin(0.01)
            if ivarSet==nVarSets-1:
                padScale=deltaY/(paddy+deltaY)
                subPads[-1].SetBottomMargin(padScale)
            else:
                subPads[-1].SetBottomMargin(0.01)
            subPads[-1].SetGridx(False)
            subPads[-1].SetGridy(True)
            subPads[-1].cd()
            garbageList.append(subPads[-1])

            ratioframe.Draw()

            allLegs.append( ROOT.TLegend(0.15,0.95,0.9,0.85) )
            allLegs[-1].SetBorderSize(0)
            allLegs[-1].SetFillStyle(0)
            allLegs[-1].SetTextFont(42)
            allLegs[-1].SetTextSize(0.1)
            if ivarSet==nVarSets-1:
                allLegs[-1].SetTextSize(0.08)
            allLegs[-1].SetNColumns(4)
            self.statUncRatioH.Draw('e2same')
            allLegs[-1].AddEntry(self.statUncRatioH,'Stat.','f')

            for i in xrange(ivarSet*2,(ivarSet+1)*2):
                if i>=len(relVars): break
                nSubVars=len(self.varH[ relVars[i].GetTitle() ])
                drawOpt='2'
                if nSubVars==1: drawOpt='p'
                relVars[i].Draw(drawOpt)
                allLegs[-1].AddEntry(relVars[i],relVars[i].GetTitle(),'f' if drawOpt!='p' else 'lp')

            allLegs[-1].Draw()

            subPads[-1].RedrawAxis()

        #all done
        p1.RedrawAxis()
        c.cd()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('%s/%s.%s'%(outdir,self.var,ext))
        for item in garbageList : item.Delete()


"""
"""
def readPlotsFrom(args,opt):

    outdir=args[0].replace('.root','')

    varTypes=OrderedDict()
    varTypes['CUET8P2MT4']=['nominal']
    for varDesc in opt.vars.split(';'):
        varTitle,varNamesEnum=varDesc.split(':')
        varList=[x for x in varNamesEnum.split(',')]
        varTypes[varTitle]=varList
    print varTypes

    fIn=ROOT.TFile.Open(args[0])
    fSyst=ROOT.TFile.Open(args[1]) if len(args)>1 else None
    for tkey in fIn.GetListOfKeys():
        key=tkey.GetName()
        if key.find('m_')==0 : continue
        if key.find('rec')<0 : continue
        if key.find('inc')<0 : continue

        #create the data, signal and total background
        t=fIn.Get(key)
        data,signal,bkg=None,None,None
        for pkey in t.GetListOfKeys():
            h=t.Get(pkey.GetName())
            if not h.InheritsFrom('TH1') : continue
            if 'Data' in h.GetTitle():
                data=h.Clone('data')
            elif h.GetTitle() in opt.signal:
                signal=h.Clone('%s nominal'%h.GetName())
            else:
                if bkg is None: bkg=h.Clone('bkg')
                else : bkg.Add(h)

        #subtract the background from the data
        try:
            data.Add(bkg,-1)
            data.Scale(1./data.Integral())
            data.SetMarkerStyle(20)
            data.SetLineColor(1)
        except:
            pass

        signal.Scale(1./signal.Integral())
        plot=VarPlot(key,varTypes)
        plot.addNominal(data)
        signal.SetTitle('CUET8P2MT4')
        plot.addVariation(signal,'nominal','CUET8P2MT4')
        for varTitle in varTypes:
            for varName in varTypes[varTitle]:
                try:
                    hvar=fSyst.Get(key).Get(signal.GetName().replace('nominal',varName))
                    print hvar,signal.GetName().replace('nominal',varName)
                    hvar.Scale(1./hvar.Integral())
                    plot.addVariation(hvar,varName,varTitle)
                except:
                    pass
        plot.show(outdir)



"""
"""
def main():

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
                      default='#deltaCUET8P2MT4:UEup,UEdn;FSR:fsr up,fsr dn;ISR:isr up,isr dn;HW++EE5C:Herwig++;m_{t}:m=166.5,m=178.5;Colour flip:cflip',
                      type='string')
    (opt, args) = parser.parse_args()


    readPlotsFrom(args,opt)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
