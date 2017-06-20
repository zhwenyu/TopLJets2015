#!/usr/bin/env/python

import ROOT
import pickle
import numpy as np
import array as array

#var name, var title, use to slice phase space, use as observable, use as event axis, is angle
VARS={
    'ptttbar'        : ('p_{T}(t#bar{t})',  True,  False, False, False),
    'ptll'           : ('p_{T}(l,l)',       True,  False, False, False),
    #'phittbar'       : ('#phi(t#bar{t})',   True,  False, True,  True),
    #    'ptpos'          : ('p_{T}(l^{+})',     True,  False, False, False),
    #    'phipos'         : ('#phi(l^{+})',      True,  False, True,  True),
    
    #'phill'          : ('#phi(ll)',         True,  False, True,  True),
    'nj'             : ('N(jets)',          True,  False, False, False),
    'chmult'         : ('N(ch)',            False,  True,  False, False),
    'chflux'         : ('#Sigma p_{T}(ch)', False, True,  False, False),
    'chavgpt'        : ('#bar{p}_{T}(ch)',  False, True,  False, False),
    'chfluxz'        : ('#Sigma p_{z}(ch)', False, True,  False, False),
    'chavgpz'        : ('#bar{p}_{z}(ch)',  False, True,  False, False),
    'sphericity'     : ('Sphericity',       False, True,  False, False),
    'aplanarity'     : ('Aplanarity',       False, True,  False, False),
    'C'              : ('C',                False, True,  False, False),
    'D'              : ('D',                False, True,  False, False)
    }

OBSVARS   = filter(lambda var : VARS[var][2], VARS)
SLICEVARS = filter(lambda var : VARS[var][1], VARS)
EVAXES    = ['ptttbar','ptll'] #filter(lambda var : VARS[var][3], VARS)

SYSTS = [ ('',   0,0,False),
          ('puup',  1,0,False),
          ('pudn',  2,0,False),
          ('effup', 3,0,False),
          ('effdn', 4,0,False),
          ('toppt', 5,0,False),
          ('murup', 9,0,False),
          ('murdn', 12,0,False),
          ('mufup', 7,0,False),
          ('mufdn', 8,0,False),
          ('qup',   10,0,False),
          ('qdn',   14,0,False),
          ('btagup',0,1,False),
          ('btagdn',0,2,False),
          ('jesup', 0,3,False),
          ('jesdn', 0,4,False),
          ('jerup', 0,5,False),
          ('jerdn', 0,6,False),
          ('eesup', 0,7,False),
          ('eesdn', 0,8,False),
          ('mesup', 0,9,False),
          ('mesdn', 0,10,False),
          ('tkeff', 0,0,1),
          ('tkeffbcdef', 0,0,2),
          ('tkeffgh', 0,0,3),
          ('tkeffeta', 0,0,4)
          ]


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
    def fillInclusive(self,obs,ue,sliceVarVals=None,ivar=0):
            
        sliceVar=None
        recSliceShift,genSliceShift=0,0
        try:
            sliceVar,genSliceVal,recSliceVal=sliceVarVals
            genSliceShift=self.getBinForVariable(genSliceVal, self.axes[ (sliceVar,False) ])-1
            recSliceShift=self.getBinForVariable(recSliceVal, self.axes[ (sliceVar,True) ])-1
        except:
            pass

        if obs==sliceVar : return
        
        #event weight
        weight=ue.w[ivar]

        #GEN level counting
        genCts=getattr(ue,'gen_chmult')
        genVal=getattr(ue,'gen_'+obs)
        genBin=self.getBinForVariable(genVal, self.axes[(obs,False)])-1
        genBin += genSliceShift*(self.axes[(obs,False)].GetNbins())            
        if not ue.gen_passSel : genBin=-1        
        if ivar==0 and genCts>0 :
            self.histos[(obs,sliceVar,'inc',None,False)].Fill(genBin,weight)

        #RECO level counting
        recCts=getattr(ue,'rec_chmult')[ivar]
        recVal=getattr(ue,'rec_'+obs)[ivar]
        recBin=self.getBinForVariable( recVal,  self.axes[(obs,True)])-1
        recBin += recSliceShift*(self.axes[(obs,True)].GetNbins())
        if not ue.rec_passSel[ivar] : recBin=-1
        if ue.rec_passSel[ivar] and recCts>0:
            if ivar==0 : 
                self.histos[(obs,sliceVar,'inc',None,True)].Fill(recBin,weight)
                if genCts==0: self.histos[(obs,sliceVar,'inc','fakes',True)].Fill(recBin,weight)
            self.histos[(obs,sliceVar,'inc','syst',True)].Fill(recBin,ivar,weight)
                
        #Migration matrix
        if genCts>0:
            key=(obs,sliceVar,'inc',ivar,'mig')
            self.histos[key].Fill(genBin,recBin,weight)

    """
    differential histogram filling
    """
    def fillDifferential(self,sliceVarVals,obs,a,ue,ivar):

        if not ue.gen_passSel and not ue.rec_passSel[ivar]: return

        sliceVar=None
        recSliceShift,genSliceShift=0,0
        try:
            sliceVar,genSliceVal,recSliceVal=sliceVarVals
            genSliceShift=self.getBinForVariable(genSliceVal, self.sliceAxes[ (sliceVar,False) ])-1
            recSliceShift=self.getBinForVariable(recSliceVal, self.sliceAxes[ (sliceVar,True) ])-1
        except:
            pass

        if obs==sliceVar : return

        #event weight   
        weight=ue.w[ivar]
        
        #GEN level counting
        genCts=getattr(ue,'gen_chmult_wrtTo')[a]        
        genVal=getattr(ue,'gen_%s_wrtTo'%obs)[a]
        genBinOffset=0
        totalBinsGen=self.axes[(obs,False)].GetNbins()
        for idx_gen in xrange(0,3):            
            if genCts[idx_gen]==0 : continue
            genBin  = idx_gen*totalBinsGen+self.getBinForVariable( genVal[idx_gen], self.axes[(obs,False)])-1
            genBin += genSliceShift*(3*totalBinsGen)
            if not ue.gen_passSel : genBin=-1
            self.histos[(obs,sliceVar,a,None,False)].Fill(genBin,weight)
            
        #RECO level counting
        recCts=getattr(ue,'rec_chmult_incWrtTo')[a]
        recVal=getattr(ue,'rec_%s_incWrtTo'%obs)[a]
        recBinOffset=0
        totalBinsRec=self.axes[(obs,True)].GetNbins()
        for idx_rec in xrange(0,3):
            if not ue.rec_passSel[ivar] : continue
            if recCts[idx_rec]==0: continue
            recBin  = idx_rec*totalBinsRec+self.getBinForVariable( recVal[idx_rec],  self.axes[(obs,True)])-1
            recBin += recSliceShift*(3*totalBinsRec)
            self.histos[(obs,sliceVar,a,None,True)].Fill(recBin,weight)
        
        #MC truth counting
#            recCtsMtrx=getattr(ue,'rec_chmult_wrtTo')[a]
#        recValMtrx=getattr(ue,'rec_%s_wrtTo'%obs)[a]
#            recBinOffset=0
#            for idx_rec in xrange(0,3):
#            
#                recBin=-1
#                if recCtsMtrx[idx_gen][idx_rec]>0:
#                    recBin  = recBinOffset+self.getBinForVariable( recVal[idx_rec],  self.axes[(obs,a,idx_rec,True)])-1
#                    recBin += recSliceShift*totalBinsRec
#                    if not passSel : recBin=-1
#                if genCts[idx_gen]>0 or recCtsMtrx[idx_gen][idx_rec]>0: 
#                    if sliceVar:
#                        self.histos[(obs,a,sliceVar)].Fill(genBin,recBin,weight)
#                    else:
#                        self.histos[(obs,a)].Fill(genBin,recBin,weight)
#
#                recBinOffset += self.axes[(obs,a,idx_rec,True)].GetNbins()
#
#            genBinOffset += self.axes[(obs,a,idx_gen,False)].GetNbins()
#


    """
    return the most appropriate bin for a given value, taking into account the range available
    """
    def getBinForVariable(self,val,axis):
        xmin,xmax=axis.GetXmin(),axis.GetXmax()
        if val>xmax : return axis.GetNbins()
        if val<xmin : return 0
        return axis.FindBin(val)
