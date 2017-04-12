#!/usr/bin/env/python

import ROOT
from TopLJets2015.TopAnalysis.eventShapeTools import *

"""
parses the event and counts particles in each region at gen/rec levels
"""
class UEEventCounter:

    """
    start variables
    """
    def __init__(self,axes=[],ptthreshold=[1.0,0.9],etathreshold=2.1,varList=[('',   0,0,False)]):

        self.ptthreshold  = ptthreshold  #pt threshold at reco levels
        self.etathreshold = etathreshold #eta threshold 
        self.piMass=0.139570

        self.varList=varList
        self.reset(axes)
        self.tkEffVar=0.5

        #tracking efficiency scale factor
        fIn=ROOT.TFile.Open('${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/era2016/MuonTracking_EfficienciesAndSF_BCDEFGH.root')
        self.tkEffSF=fIn.Get('ratio_eff_eta3_tk0_dr030e030_corr')
        fIn.Close()

            
    """
    restart counters
    """
    def reset(self,axes=None):        
        
        nvars=len(self.varList)

        self.w                  = [1.0]*nvars
        self.rec_passSel        = [False]*nvars
        self.gen_passSel        = False

        self.rec_sphericity     = [0]*nvars           #event shapes
        self.rec_aplanarity     = [0]*nvars
        self.rec_C              = [0]*nvars
        self.rec_D              = [0]*nvars
        self.gen_sphericity     = 0
        self.gen_aplanarity     = 0
        self.gen_C              = 0
        self.gen_D              = 0

        self.rec_chmult   = [0]*nvars           #inclusive rec counts
        self.rec_chflux   = [0]*nvars
        self.rec_chavgpt  = [0]*nvars
        self.rec_chfluxz  = [0]*nvars
        self.rec_chavgpz  = [0]*nvars
        self.gen_chmult   = 0                 #inclusive gen counts
        self.gen_chflux   = 0
        self.gen_chavgpt  = 0
        self.gen_chfluxz  = 0
        self.gen_chavgpz  = 0

        #self.rec_chmult_incWrtTo={}    #reco counts per reco region
        #self.rec_chflux_incWrtTo={}
        #self.rec_chavgpt_incWrtTo={}
        #self.rec_chmult_wrtTo={}       #migration matrix true-reco regions
        #self.rec_chflux_wrtTo={}
        #self.rec_chavgpt_wrtTo={}
        #self.gen_chmult_wrtTo={}       #gen counts per gen region
        #self.gen_chflux_wrtTo={}
        #self.gen_chavgpt_wrtTo={}
        #
        #if axes is None : axes=self.rec_chmult_incWrtTo.keys()
        #for a in axes:
        #    self.rec_chmult_incWrtTo[a]  = [ [0]*3 ]*nvars
        #    self.rec_chflux_incWrtTo[a]  = [ [0]*3 ]*nvars
        #    self.rec_chavgpt_incWrtTo[a] = [ [0]*3 ]*nvars
        #    self.rec_chmult_wrtTo[a]     = [[0]*3,[0]*3,[0]*3]
        #    self.rec_chflux_wrtTo[a]     = [[0]*3,[0]*3,[0]*3]
        #    self.rec_chavgpt_wrtTo[a]    = [[0]*3,[0]*3,[0]*3]
        #    self.gen_chmult_wrtTo[a]     = [0]*3
        #    self.gen_chflux_wrtTo[a]     = [0]*3
        #    self.gen_chavgpt_wrtTo[a]    = [0]*3

    """
    printout the event contents
    """
    def show(self):
        print 'Level # pTsum <pT> pZsum <pZ>'
        print 'RECO %d %3.1f %3.1f %3.1f %3.1f'%( self.rec_chmult[0], self.rec_chflux[0], self.rec_chavgpt[0], self.rec_chfluxz[0], self.rec_chavgpz[0] )
        print 'GEN %d %3.1f %3.1f %3.1f %3.1f'%( self.gen_chmult[0],  self.gen_chflux[0], self.gen_chavgpt[0], self.gen_chfluxz[0], self.gen_chavgpz[0] )

    """
    0 - tow(ards), 1 trans(verse), 2 away
    """
    def getRegionFor(self,dphi) :
        if ROOT.TMath.Abs(dphi) < ROOT.TMath.Pi()/3.     : return 0
        elif ROOT.TMath.Abs(dphi) < 2*ROOT.TMath.Pi()/3. : return 1
        return 2
    
    """
    converts index to name
    """
    def getRegionName(self,idx):
        if idx==0   : return 'tow'
        elif idx==1 : return 'tra'
        elif idx==2 : return 'awa'
        return 'inc'

    """
    count the particles in an event
    """
    def count(self,t,debug=False):
                
        self.reset()
        evshapes=EventShapeTool()
        p4=ROOT.TLorentzVector(0,0,0,0)

        for ivar in xrange(0,len(self.varList)):
            _,wgtIdx,varIdx,varyTkEff = self.varList[ivar]

            #reco level
            self.w[ivar]           = t.weight[wgtIdx]
            self.rec_passSel[ivar] = ((t.passSel>>varIdx)&0x1)
            if self.rec_passSel[ivar]:

                selP4=[]
                for n in xrange(0,t.n):

                    if t.pt[n]<self.ptthreshold[1] : continue
                    if abs(t.eta[n])>self.etathreshold : continue

                    #decide if track should still be accepted
                    if varyTkEff:
                        sf=self.tkEffSF.Eval(abs(t.eta[n]))
                        rnd=ROOT.gRandom.Uniform(1.0)
                        if rnd>(sf-(1-sf)*self.tkEffVar) : continue
                    
                        
                    p4.SetPtEtaPhiM(t.pt[n],t.eta[n],t.phi[n],self.piMass)
                    selP4.append(ROOT.TLorentzVector(p4))

                    self.rec_chmult[ivar] +=1
                    self.rec_chflux[ivar] += p4.Pt()
                    self.rec_chfluxz[ivar] += abs(p4.Pz())
                    #for a in self.rec_chmult_wrtTo:
                    # 
                    #     phirec=getattr(t,a)[varIdx]
                    #     idxrec=self.getRegionFor( ROOT.TVector2.Phi_mpi_pi(t.phi[n]-phirec) )
                    #     
                    #    phigen=getattr(t,'gen_'+a)
                    #    idxgen=self.getRegionFor( ROOT.TVector2.Phi_mpi_pi(t.phi[n]-phigen) )
                    #    self.rec_chmult_incWrtTo[a][idxrec] +=1
                    #     self.rec_chflux_incWrtTo[a][idxrec] += t.pt[n]                    
                    #     self.rec_chmult_wrtTo[a][idxgen][idxrec] +=1
                    #     self.rec_chflux_wrtTo[a][idxgen][idxrec] += t.pt[n]


                #event shapes
                evshapes.analyseNewEvent(selP4)
                self.rec_sphericity[ivar]     = evshapes.sphericity
                self.rec_aplanarity[ivar]     = evshapes.aplanarity
                self.rec_C[ivar]              = evshapes.C
                self.rec_D[ivar]              = evshapes.D
                
                #average pt/pz
                self.rec_chavgpt[ivar] = self.rec_chflux[ivar]/self.rec_chmult[ivar]  if self.rec_chmult[ivar]>0 else 0.
                self.rec_chavgpz[ivar] = self.rec_chfluxz[ivar]/self.rec_chmult[ivar] if self.rec_chmult[ivar]>0 else 0.
                #for a in self.rec_chmult_wrtTo:
                #    for k in xrange(0,len(self.rec_chmult_wrtTo[a])):
                #        for l in xrange(0,len(self.rec_chmult_wrtTo[a][k])):
                #            ncounted=self.rec_chmult_wrtTo[a][k][l]
                #            if ncounted==0 : continue
                #            self.rec_chavgpt_wrtTo[a][k][l]=self.rec_chflux_wrtTo[a][k][l]/ncounted
                #for k in xrange(0,len(self.rec_chmult_incWrtTo[a])):
                #    ncounted=self.rec_chmult_incWrtTo[a][k]
                #    if ncounted==0 : continue
                #    self.rec_chavgpt_incWrtTo[a][k]=self.rec_chflux_incWrtTo[a][k]/ncounted

        #gen level
        self.gen_passSel=(t.gen_passSel&0x1)
        if self.gen_passSel:

            selP4=[]
            for n in xrange(0,t.gen_n):

                if t.gen_pt[n] < self.ptthreshold[0] : continue
                if abs(t.gen_eta[n])>self.etathreshold : continue
                #print 'gen',n,'(%3.1f,%3.2f,%3.2f)'%(t.gen_pt[n],t.gen_eta[n],t.gen_phi[n]),t.gen_id[n]
                p4.SetPtEtaPhiM(t.gen_pt[n],t.gen_eta[n],t.gen_phi[n],self.piMass)
                selP4.append(ROOT.TLorentzVector(p4))

                self.gen_chmult +=1
                self.gen_chflux += p4.Pt()
                self.gen_chfluxz += abs(p4.Pz())
                #for a in self.gen_chmult_wrtTo:
                #                        
                #    phigen=getattr(t,'gen_'+a)
                #    idxgen=self.getRegionFor( ROOT.TVector2.Phi_mpi_pi(t.gen_phi[n]-phigen) )
                #    self.gen_chmult_wrtTo[a][idxgen] +=1
                #    self.gen_chflux_wrtTo[a][idxgen] += t.gen_pt[n]

            #event shapes
            evshapes.analyseNewEvent(selP4)
            self.gen_sphericity     = evshapes.sphericity
            self.gen_aplanarity     = evshapes.aplanarity
            self.gen_C              = evshapes.C
            self.gen_D              = evshapes.D

            #average pt/pz
            self.gen_chavgpt = self.gen_chflux/self.gen_chmult  if self.gen_chmult>0 else 0.
            self.gen_chavgpz = self.gen_chfluxz/self.gen_chmult if self.gen_chmult>0 else 0.
            #for a in self.gen_chmult_wrtTo:
            #    for k in xrange(0,len(self.gen_chmult_wrtTo[a])):
            #        ncounted=self.gen_chmult_wrtTo[a][k]
            #        if ncounted==0 : continue
            #        self.gen_chavgpt_wrtTo[a][k]=self.gen_chflux_wrtTo[a][k]/ncounted

        #print t.ptpos[0],t.sumpt[0]-t.ptpos[0]
        if debug : self.show()

