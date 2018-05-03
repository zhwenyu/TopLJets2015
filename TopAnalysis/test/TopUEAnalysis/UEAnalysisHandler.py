#!/usr/bin/env/python

import ROOT
import pickle
import numpy as np
import array as array

#var name, var title, use to slice phase space, use as observable, can be counted in regions, 
VARTITLES={
    'ptttbar'        : '|#vec{p}_{T}(t#bar{t})|',
    'ptll'           : '|#vec{p}_{T}(ll)|',
    'mll'            : 'M(l,l)',
    'nj'             : 'N(jets)',
    'chmult'         : 'N_{ch}',
    'chflux'         : '#scale[0.7]{#sum}p_{T}',
    'chavgpt'        : '#bar{p}_{T}',
    'chfluxz'        : '#scale[0.7]{#sum}p_{z}',
    'chavgpz'        : '#bar{p}_{z}',
    'chrecoil'       : '|#vec{p}_{T}|',
    'sphericity'     : 'S', #'Sphericity',
    'aplanarity'     : 'A', #'Aplanarity',
    'C'              : 'C',
    'D'              : 'D',
    'sphericity_2'     : 'S_{2}', #'Sphericity',
    'aplanarity_2'     : 'A_{2}', #'Aplanarity',
    'C_2'              : 'C_{2}',
    'D_2'              : 'D_{2}',
    'maxRap'         : 'max|#eta|',
    'rapDist'        : '#Delta#eta'
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

        #control gen-level distributions for synchronization with RIVET
        self.histos['gen_ptll']=ROOT.TH1F("gen_ptll",";Dilepton transverse momentum [GeV];Events",50,0,200)
        self.histos['gen_ptlsum']=ROOT.TH1F('gen_ptlsum',";Transverse momentum sum [GeV];Events",50,40,300) 
        self.histos['gen_mll']=ROOT.TH1F('gen_mll',";Dilepton invariant mass [GeV];Events",50,0,400)
        self.histos['gen_ptpos']=ROOT.TH1F("gen_ptpos",";l^{+} transverse momentum [GeV];Events",50,0,200)
        self.histos['gen_nj']=ROOT.TH1F('gen_nj',";Extra jet multiplicity;Events",7,0,7)
        self.histos['gen_bpt']=ROOT.TH1F('gen_bpt',";b p_{T};Events",50,0,250)        
        self.histos['gen_jpt']=ROOT.TH1F('gen_jpt',";Extra jet p_{T};Events",50,0,250)
  
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
        if not ue.gen_passSel : 
            genBin=-1
        else:
            self.histos['gen_ptll'].Fill(ue.gen_ptll)
            self.histos['gen_ptlsum'].Fill(ue.gen_sumpt)
            self.histos['gen_mll'].Fill(ue.gen_mll)
            self.histos['gen_ptpos'].Fill(ue.gen_ptpos)
            self.histos['gen_nj'].Fill(ue.gen_nj)
            for pt in ue.gen_jpt:
                if pt==0 : continue
                self.histos['gen_jpt'].Fill(pt)
            for pt in ue.gen_bpt:
                self.histos['gen_bpt'].Fill(pt)            
            if genCts>0:
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
            if ue.gen_passSel and genCts>0:
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
