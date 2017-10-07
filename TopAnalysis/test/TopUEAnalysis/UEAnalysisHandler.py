#!/usr/bin/env/python

import ROOT
import pickle
import numpy as np
import array as array

#var name, var title, use to slice phase space, use as observable, can be counted in regions, 
VARTITLES={
    'ptttbar'        : 'p_{T}(t#bar{t})',
    'ptll'           : 'p_{T}(l,l)',
    'nj'             : 'N(jets)',
    'chmult'         : 'N(ch)',
    'chflux'         : '#Sigma p_{T}(ch)',
    'chavgpt'        : '#bar{p}_{T}(ch)',
    'chfluxz'        : '#Sigma p_{z}(ch)',
    'chavgpz'        : '#bar{p}_{z}(ch)',
    'sphericity'     : 'Sphericity',
    'aplanarity'     : 'Aplanarity',
    'C'              : 'C',
    'D'              : 'D'
    }

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
          ('tkeffeta', 0,0,4),
          ('tkeffdstar', 0,0,5)
          ]


"""
parses the event and counts particles in each region at gen/rec levels
"""
class UEAnalysisHandler:

    def __init__(self,analysisCfg,useSysts):
        
        with open(analysisCfg,'r') as cachefile:
            self.analysisCfg = pickle.load(cachefile)
            self.cuts        = pickle.load(cachefile)
            self.obs         = pickle.load(cachefile)
            self.ptthreshold = pickle.load(cachefile)

        self.nVars=len(SYSTS) if useSysts else 1
        self.histos={}
        self.histos['gen']=self.analysisCfg[('gen','histo')].Clone('gen')
        for i in xrange(0,self.nVars):
            self.histos['reco_%d'%i]=self.analysisCfg[('reco','histo')].Clone('reco_%d'%i)
            self.histos['reco_%d'%i].SetTitle(SYSTS[i][0])
            self.histos['fakes_%d'%i]=self.analysisCfg[('reco','histo')].Clone('fakes_%d'%i)
            self.histos['fakes_%d'%i].SetTitle(SYSTS[i][0])
            self.histos['mig_%d'%i]=self.analysisCfg[('mig','histo')].Clone('mig_%d'%i)
            self.histos['mig_%d'%i].SetTitle(SYSTS[i][0])
        for key in self.histos:
            self.histos[key].Sumw2()
            self.histos[key].SetDirectory(0)

        #print out
        print '[UEAnalysisHandler] configured with'
        print '\t obs=',self.obs
        print '\t ',len(self.histos),'histos to fill, corresponding to',self.nVars,'variations'
        print '\t cuts=',self.cuts
        print '\t kinematics=',self.ptthreshold

    """
    inclusive histogram filling
    """
    def fillHistos(self,ue):
        
        #GEN level counting
        genCts=getattr(ue,'gen_chmult')
        genVal=getattr(ue,'gen_'+self.obs)  
        genBin=self.getBinForVariable(genVal, self.analysisCfg[('gen','axis')])-1
        if not ue.gen_passSel : genBin=-1               
        if genCts>0 :
            self.histos['gen'].Fill(genBin,ue.w[0])

        #RECO level counting (loop over variations)
        for i in xrange(0,self.nVars):
        
            #event weight
            weight=ue.w[i]
        
            #RECO level counting
            recCts=getattr(ue,'rec_chmult')[i]
            recVal=getattr(ue,'rec_'+self.obs)[i]
            recBin=self.getBinForVariable( recVal, self.analysisCfg[('reco','axis')])-1
            if not ue.rec_passSel[i] : recBin=-1            
            if ue.rec_passSel[i] and recCts>0:
                self.histos['reco_%d'%i].Fill(recBin,weight)
                if not ue.gen_passSel : 
                    self.histos['fakes_%d'%i].Fill(recBin,weight)
 
            #Migration matrix
            if genCts>0:
                self.histos['mig_%d'%i].Fill(genBin,recBin,weight)


    """
    return the most appropriate bin for a given value, taking into account the range available
    """
    def getBinForVariable(self,val,axis):
        xmin,xmax=axis.GetXmin(),axis.GetXmax()       
        if val>=xmax : return axis.GetNbins()
        if val<xmin : return 0
        xbin=axis.FindBin(val)
        return xbin
