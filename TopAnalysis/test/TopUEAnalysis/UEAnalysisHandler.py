#!/usr/bin/env/python

import ROOT
import pickle
import numpy as np
import array as array

#var name, var title, use to slice phase space, use as observable, can be counted in regions, 
VARS={
    'ptttbar'        : ('p_{T}(t#bar{t})',  True,  False, False),
    'ptll'           : ('p_{T}(l,l)',       True,  False, False),
    'nj'             : ('N(jets)',          True,  False, False),
    'chmult'         : ('N(ch)',            False, True,  True),
    'chflux'         : ('#Sigma p_{T}(ch)', False, True,  True),
    'chavgpt'        : ('#bar{p}_{T}(ch)',  False, True,  True),
    'chfluxz'        : ('#Sigma p_{z}(ch)', False, True,  False),
    'chavgpz'        : ('#bar{p}_{z}(ch)',  False, True,  False),
    'sphericity'     : ('Sphericity',       False, True,  False),
    'aplanarity'     : ('Aplanarity',       False, True,  False),
    'C'              : ('C',                False, True,  False),
    'D'              : ('D',                False, True,  False),
    }

OBSVARS   = filter(lambda var : VARS[var][2], VARS)
SLICEVARS = filter(lambda var : VARS[var][1], VARS)
EVAXES    = ['ptttbar','ptll']

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
    def fillHistos(self,obs,ue,sliceVarVals=None,ivar=0):
            
        sliceVar=None
        recSliceBins,genSliceBins=[0],[0]
        axes=['inc']
        try:
            sliceVar,genSliceVal,recSliceVal=sliceVarVals
            recSliceBins.append( self.getBinForVariable(recSliceVal, self.axes[ (sliceVar,False) ]) )
            genSliceBins.append( self.getBinForVariable(genSliceVal, self.axes[ (sliceVar,False) ]) )
            if sliceVar in EVAXES and VARS[obs][3] : axes+=[0,1,2]
        except:
            pass

        if obs==sliceVar : return
        
        #event weight
        weight=ue.w[ivar]
        
        for axis in axes:

            #GEN level counting
            genCts=getattr(ue,'gen_chmult') if axis=='inc' else getattr(ue,'gen_chmult_wrtTo')[sliceVar][axis]
            genVal=getattr(ue,'gen_'+obs)   if axis=='inc' else getattr(ue,'gen_%s_wrtTo'%obs)[sliceVar][axis]
            genBin=self.getBinForVariable(genVal, self.axes[(obs,False)])-1
            if not ue.gen_passSel : genBin=-1        
            if ivar==0 and genCts>0 :
                for b in genSliceBins:
                    self.histos[(obs,sliceVar,b,axis,None,False)].Fill(genBin,weight)

            #RECO level counting
            recCts=getattr(ue,'rec_chmult')[ivar] if axis=='inc' else getattr(ue,'rec_chmult_wrtTo')[sliceVar][ivar][axis]
            recVal=getattr(ue,'rec_'+obs)[ivar]   if axis=='inc' else getattr(ue,'rec_%s_wrtTo'%obs)[sliceVar][ivar][axis]
            recBin=self.getBinForVariable( recVal,  self.axes[(obs,True)])-1
            if not ue.rec_passSel[ivar] : recBin=-1
            if ue.rec_passSel[ivar] and recCts>0:
                for b in recSliceBins:
                    if ivar==0 : 
                        self.histos[(obs,sliceVar,b,axis,None,True)].Fill(recBin,weight)
                        if genCts==0: self.histos[(obs,sliceVar,b,axis,'fakes',True)].Fill(recBin,weight)
                    self.histos[(obs,sliceVar,b,'inc','syst',True)].Fill(recBin,ivar,weight)
                
            #Migration matrix
            if genCts>0:
                for b in genSliceBins:
                    key=(obs,sliceVar,b,axis,ivar,'mig')
                    self.histos[key].Fill(genBin,recBin,weight)


    """
    return the most appropriate bin for a given value, taking into account the range available
    """
    def getBinForVariable(self,val,axis):
        xmin,xmax=axis.GetXmin(),axis.GetXmax()       
        if val>=xmax : return axis.GetNbins()
        if val<xmin : return 0
        return axis.FindBin(val)
