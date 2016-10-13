#!/usr/bin/env/python

import ROOT
import pickle
import numpy as np
import array as array


"""
parses the event and counts particles in each region at gen/rec levels
"""
class UEAnalysisHandler:

    def __init__(self,analysisCfg):
        
        with open(analysisCfg,'r') as cachefile:
            self.axes=pickle.load(cachefile)
            self.histos=pickle.load(cachefile)

    """
    inclusive histogram filling
    """
    def fillInclusive(self,obs,ue,weight,gen_passSel,passSel):

        if not passSel and not gen_passSel: return

        genCts=getattr(ue,'gen_chmult')
        genVal=getattr(ue,'gen_'+obs)
        genBin=self.getBinForVariable(genVal, self.axes[(obs,'inc',0,False)])-1
        if gen_passSel and genCts>0: 
            self.histos[(obs,'inc',False)].Fill(genBin,weight)
                
        recCts=getattr(ue,'rec_chmult')
        recVal=getattr(ue,'rec_'+obs)
        recBin=self.getBinForVariable( recVal,  self.axes[(obs,'inc',0,True)])-1
        if passSel and recCts>0:
            self.histos[(obs,'inc',True)].Fill(recBin,weight)
        
        #underflows will be used for events failing gen/reco selections
        if not gen_passSel : genBin=-1
        if not passSel     : recBin=-1
        if recCts>0 : self.histos[(obs,'inc')].Fill(genBin,recBin,weight)

    """
    differential histogram filling
    """
    def fillDifferential(self,obs,a,ue,weight,gen_passSel,passSel):

        if not passSel and not gen_passSel: return

        recCtsMtrx=getattr(ue,'rec_chmult_wrtTo')[a]
        genCtsMtrx=getattr(ue,'gen_chmult_wrtTo')[a]

        recMtrx=getattr(ue,'rec_%s_wrtTo'%obs)[a]
        genMtrx=getattr(ue,'gen_%s_wrtTo'%obs)[a]

        genBinOffset=0
        for idx_gen in xrange(0,3):

            genCts=genCtsMtrx[idx_gen]
            genVal=genMtrx[idx_gen]
            genBin=genBinOffset+self.getBinForVariable( genVal, self.axes[(obs,a,idx_gen,False)])-1
            if gen_passSel and genCts>0 : 
                self.histos[(obs,a,False)].Fill(genBin,weight)
                
            recBinOffset=0
            for idx_rec in xrange(0,3):

                recCts=recCtsMtrx[idx_gen][idx_rec]
                recVal=recMtrx[idx_gen][idx_rec]            
                recBin=recBinOffset+self.getBinForVariable( recVal,  self.axes[(obs,a,idx_rec,True)])-1
                
                if passSel and recCts>0 : 
                    self.histos[(obs,a,True)].Fill(recBin,weight)

                if not passSel     : recBin=-1
                if not gen_passSel : genBin=-1
                if recCts>0        : self.histos[(obs,a)].Fill(genBin,recBin,weight)

                recBinOffset += self.axes[(obs,a,idx_rec,True)].GetNbins()

            genBinOffset += self.axes[(obs,a,idx_gen,False)].GetNbins()


    """
    return the most appropriate bin for a given value, taking into account the range available
    """
    def getBinForVariable(self,val,axis):
        xmin,xmax=axis.GetXmin(),axis.GetXmax()
        if val>xmax : return axis.GetNbins()
        if val<xmin : return 0
        return axis.FindBin(val)
