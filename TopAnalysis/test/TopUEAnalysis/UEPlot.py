#!/usr/bin/env python

import ROOT
import numpy as np
from collections import defaultdict

class UEPlot:
    """A container for a graph with variations. A plot has the following properties:
    name, title, central mean and associated uncertainties, central graph and associated uncertainties,
    list of variation histograms and associated list of means
    """

    def __init__(self,name,title,ci,marker,fill,trueAxis):
        """Setup the attributes of this plot"""
        self.name=name
        self.mean=[0.,[]]
        self.plot=[ROOT.TGraphErrors(),[]]
        self.plot[0].SetName('%s_central'%self.name)
        self.plot[0].SetTitle(title)
        self.plot[0].SetMarkerStyle(marker)
        self.plot[0].SetMarkerColor(ci)
        self.plot[0].SetLineColor(ci)
        self.plot[0].SetFillColor(ci)
        self.plot[0].SetFillStyle(fill)
        for i in ['stat','exp','th']:
            self.plot[1].append( self.plot[0].Clone( '%s_%s'%(self.name,i) ) )
            self.mean[1].append( 0. )
        self.trueAxis=trueAxis        
        self.variations=defaultdict(list)
        self.variationMeans=defaultdict(list)

    def addVariation(self,varName,varType,varH):
        """Add one more variation. varType is used to signal if is nominal (None), experimental (exp) or theory (th)"""

        if varH.Integral()<=0 : return

        key=(varName,varType)
        cloneName='%s_%s_%s_%d'%(self.name,varName,varType,len(self.variations[key]))
        self.variations[key].append(varH.Clone(cloneName))
        self.variations[key][-1].SetDirectory(0)

        x,w=[],[]
        for xbin in xrange(1,self.trueAxis.GetNbins()+1):
            x.append( self.trueAxis.GetBinCenter(xbin) )
            w.append( varH.GetBinContent(xbin) )
        try:
            avg=np.average(x,weights=w)
            dx2=[(k-avg)**2 for k in x]
            var = np.average(dx2, weights=w)
            self.variationMeans[key].append( (avg,ROOT.TMath.Sqrt(var/sum(w))) )
        except:
            self.variationMeans[key].append( None )
            

    def clear(self):
        """Free memory of this plot"""
        self.plot[0].Delete()
        for p in self.plot[1]:
            p.Delete()
        for key in self.variations:
            for p in self.variations[key]:
                p.Delete()

    def finalize(self):
        """To be used once all variations have been given: fills plots and computes means"""

        #start by organizing the keys
        nomKey=None
        expKeys=[]
        thKeys=[]
        for key in self.variations:
            if key[1] is None : nomKey=key
            elif key[1]=='exp': expKeys.append(key)
            elif key[1]=='th':  thKeys.append(key)

        #check normalization
        norm=self.variations[nomKey][0].Integral()
        if norm==0: raise Exception('Null counts found in %s plot - unable to finalize'%self.name)

        for xbin in xrange(1,self.trueAxis.GetNbins()+1):
            cen        = self.trueAxis.GetBinCenter(xbin)
            xwid       = self.trueAxis.GetBinWidth(xbin)*0.5
            cts        = self.variations[nomKey][0].GetBinContent(xbin)
            
            cenVal=cts/(xwid*norm)
            self.plot[0].SetPoint(xbin-1, cen, cenVal)

            #statistical uncertainty
            ctsStatUnc = self.variations[nomKey][0].GetBinError(xbin)
            self.plot[1][0].SetPoint(xbin-1, cen, cenVal)
            self.plot[1][0].SetPointError(xbin-1, xwid, ctsStatUnc/(xwid*norm) )

            #experimental uncertainty
            ctsExpUnc=0
            for key in expKeys:
                iCtsExpUnc=0
                for h in self.variations[key]:
                    iCtsExpUnc=max( abs(h.GetBinContent(xbin)-cts), iCtsExpUnc )
                ctsExpUnc += iCtsExpUnc**2
            ctsExpUnc=ROOT.TMath.Sqrt(ctsExpUnc)
            self.plot[1][1].SetPoint(xbin-1, cen, cenVal)
            self.plot[1][1].SetPointError(xbin-1, xwid, ctsExpUnc/(xwid*norm) )
            
            #theory uncertainty
            ctsThUnc=0
            for key in thKeys:
                iCtsThUnc=0
                for h in self.variations[key]:
                    iCtsThUnc=max( abs(h.GetBinContent(xbin)-cts), iCtsThUnc )
                ctsThUnc += iCtsExpUnc**2
            ctsThUnc=ROOT.TMath.Sqrt(ctsThUnc)
            self.plot[1][2].SetPoint(xbin-1, cen, cenVal)
            self.plot[1][2].SetPointError(xbin-1, xwid, ctsThUnc/(xwid*norm) )
            
            ctsTotalUnc=ROOT.TMath.Sqrt(ctsStatUnc**2+ctsExpUnc**2+ctsThUnc**2)
            self.plot[0].SetPointError(xbin,xwid,ctsTotalUnc)

        #final mean
        self.mean[0]=self.variationMeans[nomKey][0][0]
        self.mean[1][0]=self.variationMeans[nomKey][0][1]

        #experimental uncertainty
        meanExpUnc=0
        for key in expKeys:
            iMeanExpUnc=0
            for m,_ in self.variationMeans[key]:
                iMeanExpUnc=max( abs(m-self.mean[0]), iMeanExpUnc )
            meanExpUnc += iMeanExpUnc**2
        self.mean[1][1]=ROOT.TMath.Sqrt(meanExpUnc)

        #theory uncertainty
        meanThUnc=0
        for key in thKeys:
            iMeanThUnc=0
            for m,_ in self.variationMeans[key]:
                iMeanThUnc=max( abs(m-self.mean[0]), iMeanThUnc )
            meanThUnc += iMeanThUnc**2
        self.mean[1][2]=ROOT.TMath.Sqrt(meanThUnc)

