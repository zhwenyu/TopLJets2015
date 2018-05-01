#!/usr/bin/env python

import ROOT
import numpy as np
from collections import defaultdict
from UETools import formatGraph

class UEPlot:
    """A container for a graph with variations. A plot has the following properties:
    name, title, central mean and associated uncertainties, central graph and associated uncertainties,
    list of variation histograms and associated list of means
    """

    def __init__(self,name,title,trueAxis):
        """Setup the attributes of this plot"""
        self.name=name
        self.mean=[0.,[]]
        self.meanUncTable=''
        self.plot=[ROOT.TGraphErrors(),[]]
        self.plot[0].SetName('%s_central'%self.name)
        self.plot[0].SetTitle(title)
        for i in ['stat','exp','th']:
            self.plot[1].append( self.plot[0].Clone( '%s_%s'%(self.name,i) ) )
            self.mean[1].append( 0. )
        self.trueAxis=trueAxis        
        self.variations=defaultdict(list)
        self.normVals=defaultdict(list)
        self.variationMeans=defaultdict(list)
        self.covMatrices=defaultdict(list)
        self.relUncertaintyH={}        

    def addVariation(self,varName,varType,varH):
        """Add one more variation. varType is used to signal if is nominal (None), experimental (exp) or theory (th)"""
        
        key=(varName,varType)
        try:
            if varH.Integral()<=0 : return
        except:
            print "Unable to evaluate integral for ",key
            return
        
        cloneName='%s_%s_%s_%d'%(self.name,varName,varType,len(self.variations[key]))
        self.variations[key].append(varH.Clone(cloneName))
        self.variations[key][-1].SetDirectory(0)

        #
        #compute mean of this distribution
        #
        #uncorrelated case
        x,w,w2=[],[],[]
        for xbin in xrange(1,self.trueAxis.GetNbins()+1):
            x.append( self.trueAxis.GetBinCenter(xbin) )
            w.append( varH.GetBinContent(xbin) )
            w2.append( varH.GetBinContent(xbin)**2 )
        try:
            avg=np.average(x,weights=w)
            dx2=[(k-avg)**2 for k in x]
            var = np.average(dx2, weights=w)
            neff=sum(w)
            #neff=(sum(w)**2)/sum(w2)
            self.variationMeans[key].append( (avg,ROOT.TMath.Sqrt(var/neff)) )
                
            #avg2=np.average(x**2,weights=w)
            #dx4=[(k**2-avg)**2 for k in x]
            #var = np.average(dx4, weights=w)
        except:
            self.variationMeans[key].append( None )

        #start a relative uncertainty histogram for a new variation
        if varName in self.relUncertaintyH: return
        self.relUncertaintyH[varName]=varH.Clone('%s_relunc'%cloneName)
        self.relUncertaintyH[varName].SetDirectory(0)
        self.relUncertaintyH[varName].Reset('ICE')
        
        
    def clear(self):
        """Free memory of this plot"""
        self.plot[0].Delete()
        for p in self.plot[1]:
            p.Delete()
        for key in self.variations:
            for p in self.variations[key]:
                p.Delete()

    def finalizeCovMatrices(self,statCov=None):
        """finalize covariance matrices: if received, statCovariance is assumed not-normalized"""
        
        #find the nominal prediction
        nomKey=None
        for key in self.variations:
            if not key[1] is None : continue 
            nomKey=key
            break
        if nomKey is None: return
        nbins = self.variations[nomKey][0].GetNbinsX()
        norm  = self.variations[nomKey][0].Integral()

        #systematic covariance
        #see https://twiki.cern.ch/twiki/bin/view/CMS/TopUnfolding#Treatment_of_systematic_uncertai
        totalSyst=ROOT.TMatrixF(nbins,nbins)
        for key in self.variations:
            if key ==nomKey : continue
            if len(self.variations[key])==0 : continue

            varNorms=[x.Integral() for x in self.variations[key]]
            cov=ROOT.TMatrixF(nbins,nbins)

            #build the arrays of differences between the variations and the nominal result
            dyLists=[]
            for i in xrange(0,len(self.variations[key])):
                dyLists.append(
                    [ self.variations[key][i].GetBinContent(xbin)/varNorms[i]-self.variations[nomKey][0].GetBinContent(xbin)/norm for xbin in xrange(1,nbins+1) ]
                    )

            #loop to maximize the difference at each bin
            for i in xrange(1,nbins+1):

                for j in xrange(1,nbins+1):

                    #single sided
                    maxdx,maxdy,sign=dyLists[0][i-1],dyLists[0][j-1],1.0
                    
                    #double sided
                    if len(dyLists)==2:
                        maxdx = max(abs(dyLists[0][i-1]),abs(dyLists[1][i-1]))
                        maxdy = max(abs(dyLists[0][j-1]),abs(dyLists[1][j-1]))
                        sign  = np.sign(dyLists[0][i-1]-dyLists[1][i-1])*np.sign(dyLists[0][j-1]-dyLists[1][j-1])
                        
                    #>2: take the envelope and proceed as single sided
                    if len(dyLists)>2:
                        maxdx,maxdy=0,0
                        for k in xrange(0,len(dyLists)):
                            maxdx=max(abs(dyLists[k][i-1]),maxdx)
                            maxdy=max(abs(dyLists[k][j-1]),maxdy)

                    cov[i-1][j-1]=maxdx*maxdy*sign

            self.covMatrices[key[0]]=cov
            totalSyst+=cov
        self.covMatrices['syst']=totalSyst

        #build the total covariance matrix
        self.covMatrices['total']=None
        for key in self.covMatrices:
            if key=='total' : continue
            if self.covMatrices['total'] is None:
                self.covMatrices['total']=ROOT.TMatrixF(self.covMatrices[key])
            else:
                self.covMatrices['total']+=self.covMatrices[key]

        #show for debug
        #ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
        #ROOT.gStyle.SetPaintTextFormat(".1g");
        #c=ROOT.TCanvas('c','c',500,500)
        #c.SetRightMargin(0.2)
        #c.Print('~/www/cov.pdf[')
        #for key in self.covMatrices:
        #    print key
        #    self.covMatrices[key].Draw('colztext')            
        #    c.BuildLegend(0.6,0.9,0.9,0.8,key)
        #    c.Print('~/www/cov.pdf')
        #inv=self.covMatrices['total'].InvertFast()
        #inv.Draw('colztext')
        #c.BuildLegend(0.6,0.9,0.9,0.8,'inverse')
        #c.Print('~/www/cov.pdf')
        #c.Print('~/www/cov.pdf]')
        #raise


    def finalize(self,doCov=False):
        """To be used once all variations have been given: fills plots and computes means"""

        #start by organizing the keys and normalizing distributions
        nomKey=None
        expKeys=[]
        thKeys=[]        
        for key in self.variations:
            if key[1] is None : nomKey=key
            elif key[1]=='exp': expKeys.append(key)
            elif key[1]=='th':  thKeys.append(key)
            for h in self.variations[key]: 
                hinteg=h.Integral()
                self.normVals[key].append(hinteg)
                h.Scale(1./hinteg)

        #check normalization
        if self.normVals[nomKey][0]==0: raise Exception('Null counts found in %s plot - unable to finalize'%self.name)
        for xbin in xrange(1,self.trueAxis.GetNbins()+1):
            cen        = self.trueAxis.GetBinCenter(xbin)
            xwid       = self.trueAxis.GetBinWidth(xbin)*0.5
            cts        = self.variations[nomKey][0].GetBinContent(xbin)
            
            cenVal=cts/xwid
            self.plot[0].SetPoint(xbin-1, cen, cenVal)

            #statistical uncertainty
            ctsStatUnc = self.variations[nomKey][0].GetBinError(xbin)
            self.plot[1][0].SetPoint(xbin-1, cen, cenVal)
            self.plot[1][0].SetPointError(xbin-1, xwid, ctsStatUnc/xwid )
            self.relUncertaintyH[nomKey[0]].SetBinContent(xbin,ctsStatUnc/cts if cenVal>0 else 0)


            #experimental uncertainty
            ctsExpUnc=0
            for key in expKeys:
                iCtsExpUnc=0
                for h in self.variations[key]:
                    iCtsExpUnc=max( abs(h.GetBinContent(xbin)-cts), iCtsExpUnc )
                ctsExpUnc += iCtsExpUnc**2
                self.relUncertaintyH[key[0]].SetBinContent(xbin,iCtsExpUnc/cts if cenVal>0 else 0)
            ctsExpUnc=ROOT.TMath.Sqrt(ctsExpUnc)
            self.plot[1][1].SetPoint(xbin-1, cen, cenVal)
            self.plot[1][1].SetPointError(xbin-1, xwid, ctsExpUnc/xwid )
            
            #theory uncertainty
            ctsThUnc=0
            for key in thKeys:
                iCtsThUnc=0
                for h in self.variations[key]:
                    iCtsThUnc=max( abs(h.GetBinContent(xbin)-cts), iCtsThUnc )
                ctsThUnc += iCtsThUnc**2
                self.relUncertaintyH[key[0]].SetBinContent(xbin,iCtsThUnc/cts if cenVal>0 else 0)
            ctsThUnc=ROOT.TMath.Sqrt(ctsThUnc)
            self.plot[1][2].SetPoint(xbin-1, cen, cenVal)
            self.plot[1][2].SetPointError(xbin-1, xwid, ctsThUnc/xwid)

            ctsTotalUnc=ROOT.TMath.Sqrt(ctsStatUnc**2+ctsExpUnc**2+ctsThUnc**2)
            self.plot[0].SetPointError(xbin-1,xwid,ctsTotalUnc/xwid)

        #compute covariance matrcies
        if doCov: self.finalizeCovMatrices()

        #final mean
        self.mean[0]=self.variationMeans[nomKey][0][0]
        self.mean[1][0]=self.variationMeans[nomKey][0][1]

        self.meanUncTable='Value & %3.3f\nStat & %3.2f\n'%(self.mean[0],self.mean[1][0]/self.mean[0]*100)

        #experimental uncertainty
        meanExpUnc=0
        for key in expKeys:
            iMeanExpUnc=0
            for m,_ in self.variationMeans[key]:
                iMeanExpUnc=max( abs(m-self.mean[0]), iMeanExpUnc )
            meanExpUnc += iMeanExpUnc**2
            self.meanUncTable+='%s & %3.2f\n'%(key[0],iMeanExpUnc/self.mean[0]*100)
        self.mean[1][1]=ROOT.TMath.Sqrt(meanExpUnc)

        #theory uncertainty
        meanThUnc=0
        for key in thKeys:
            iMeanThUnc=0
            for m,_ in self.variationMeans[key]:
                iMeanThUnc=max( abs(m-self.mean[0]), iMeanThUnc )
            meanThUnc += iMeanThUnc**2
            self.meanUncTable+='%s & %3.2f\n'%(key[0],iMeanThUnc/self.mean[0]*100)
        self.mean[1][2]=ROOT.TMath.Sqrt(meanThUnc)

        totalUnc=ROOT.TMath.Sqrt(meanThUnc+meanExpUnc+self.mean[1][0]**2)
        self.meanUncTable+='Total & %3.2f'%(totalUnc/self.mean[0]*100)

    def format(self,fill,color,marker,keepXUnc,shiftX):
        """
        Apply a common format to the plots
        """
        for p in self.plot[1]+[self.plot[0]]: formatGraph(p,fill,color,marker,keepXUnc,shiftX)

def getRatiosWithRespectTo(uePlots,refKey,uncType=''):
    """
    Compute the ratio using a given key as reference
    """
    uePlotRatios={}
    x,y=ROOT.Double(0),ROOT.Double(0)
    xref,yref=ROOT.Double(0),ROOT.Double(0)    
    for key in uePlots:
        numPlot=uePlots[key].plot[0]
        denPlot=uePlots[refKey].plot[0]
        if uncType=='stat' : 
            numPlot=uePlots[key].plot[1][0]
            denPlot=uePlots[refKey].plot[1][0]

        uePlotRatios[key]=numPlot.Clone('%s_2_%s_ratio%s'%(uePlots[key].name,uePlots[refKey].name,uncType))

        for np in xrange(0,numPlot.GetN()):
            
            numPlot.GetPoint(np,x,y)
            denPlot.GetPoint(np,xref,yref)
            ratio=-1 if yref==0 else y/yref

            ex=numPlot.GetErrorX(np)
            ey=numPlot.GetErrorY(np)
            #eyref=denPlot.GetErrorY(np)
                
            ratioUnc=(ey*yref)**2
            #if key!=refKey:ratioUnc+=(eyref*y)**2
            ratioUnc=0 if yref==0 else ROOT.TMath.Sqrt(ratioUnc)/(yref**2)

            uePlotRatios[key].SetPoint(np,x,ratio)
            uePlotRatios[key].SetPointError(np,ex,ratioUnc)

    return uePlotRatios
