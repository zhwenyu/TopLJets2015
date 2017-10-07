#!/usr/bin/env/python

import ROOT
from TopLJets2015.TopAnalysis.eventShapeTools import *

AXISANGLE = {
    'ptttbar':'phittbar',
    'ptll':'phill',
    'ptpos':'phipos'
    }

class UEEventCounter:
    """
    parses the event and counts particles in each region at gen/rec levels
    """
    def __init__(self,
                 ptthreshold=[0.9,0.9],
                 etathreshold=2.1,
                 cuts={},
                 systList=[('',   0,0,False)]):
        """
        start class
        """

        self.ptthreshold  = ptthreshold  #pt threshold 
        self.etathreshold = etathreshold #eta threshold 
        self.piMass       = 0.139570     #charged pion mass
        self.cuts         = cuts         #further cuts
        self.systList     = systList     #list of systematic variations

        #tracking efficiency scale factors
        self.tkEffSF={}
        for era in ['BCDEF','GH']:

            #muon based scale factors
            fIn=ROOT.TFile.Open('${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/era2016/MuonTracking_EfficienciesAndSF_%s.root'%era)
            self.tkEffSF[era]={'abseta':fIn.Get('ratio_eff_aeta_tk0_dr030e030_corr'),
                               'vtx':fIn.Get('ratio_eff_vtx_tk0_dr030e030_corr')}
            fIn.Close()

            #dstar scale factors
            self.tkEffSF[era+"_dstar"]={'eta':ROOT.TGraphErrors()}
            dstarSF=[(1.95,0.45,1.12,0.05),(-1.15,0.35,1.07,0.07),(0,0.8,1.04,0.03),(1.15,0.35,1.07,0.06),(1.95,0.45,1.12,0.05)]
            if era=='BDCDEF':
                dstarSF=[(1.95,0.45,0.93,0.04),(-1.15,0.35,1.08,0.04),(0,0.8,1.01,0.03),(1.15,0.35,1.08,0.04),(1.95,0.45,0.93,0.04)]
            for eta,etaunc,sf,sfunc in dstarSF :
                np=self.tkEffSF[era+"_dstar"]['eta'].GetN()
                self.tkEffSF[era+"_dstar"]['eta'].SetPoint(np,eta,sf)
                self.tkEffSF[era+"_dstar"]['eta'].SetPointError(np,etaunc,sfunc)


        #aux. class to compute event shapes
        self.evshapes=EventShapeTool()

        #reset counters
        self.reset()

    def reset(self):
        """
        reset counters
        """
    
        nSysts=len(self.systList)

        self.w                  = [1.0]*nSysts    #event weights
        self.rec_passSel        = [False]*nSysts  #event selection
        self.gen_passSel        = False

        self.rec_sphericity     = [0]*nSysts      #event shapes
        self.rec_aplanarity     = [0]*nSysts
        self.rec_C              = [0]*nSysts
        self.rec_D              = [0]*nSysts
        self.gen_sphericity     = 0
        self.gen_aplanarity     = 0
        self.gen_C              = 0
        self.gen_D              = 0

        self.rec_chmult   = [0]*nSysts            #classic observables
        self.rec_chflux   = [0]*nSysts
        self.rec_chavgpt  = [0]*nSysts
        self.rec_chfluxz  = [0]*nSysts
        self.rec_chavgpz  = [0]*nSysts
        self.gen_chmult   = 0 
        self.gen_chflux   = 0
        self.gen_chavgpt  = 0
        self.gen_chfluxz  = 0
        self.gen_chavgpz  = 0


    def show(self):
        """
        printout the event contents
        """
        print 'Level # pTsum <pT> pZsum <pZ>'
        print 'RECO %d %3.1f %3.1f %3.1f %3.1f'%( self.rec_chmult[0], self.rec_chflux[0], self.rec_chavgpt[0], self.rec_chfluxz[0], self.rec_chavgpz[0] )
        print 'GEN %d %3.1f %3.1f %3.1f %3.1f'%( self.gen_chmult[0],  self.gen_chflux[0], self.gen_chavgpt[0], self.gen_chfluxz[0], self.gen_chavgpz[0] )
        
    def getRegionFor(self,dphi) :
        """
        0 - tow(ards), 1 trans(verse), 2 away
        """
        if ROOT.TMath.Abs(dphi) < ROOT.TMath.Pi()/3.     : return 0
        elif ROOT.TMath.Abs(dphi) < 2*ROOT.TMath.Pi()/3. : return 1
        return 2
    
    def getRegionName(self,idx):
        """
        converts index to name
        """
        if idx==0   : return 'tow'
        elif idx==1 : return 'tra'
        elif idx==2 : return 'awa'
        return 'inc'

    def count(self,t,isMC=False,debug=False):
        """
        count the particles in an event
        """
        
        #assign an era randomly (only used for MC)
        mceraLumi=ROOT.gRandom.Uniform(35874.8)
        mcera='BCDEF' if mceraLumi < 19323.4 else 'GH'

        self.reset()

        p4=ROOT.TLorentzVector(0,0,0,0)

        #repeat N-selections for every possible variations
        for iSystVar in xrange(0,len(self.systList)):

            #weight, variation index in tree and tk. eff flag
            _,wgtIdx,varIdx,varyTkEff = self.systList[iSystVar]


            #reco level
            self.w[iSystVar]           = t.weight[wgtIdx]
            self.rec_passSel[iSystVar] = ((t.passSel>>varIdx)&0x1)

            passExtraCuts=True
            for cutKey in self.cuts:
                if cutKey=='region': continue
                cutVar=cutKey
                cutVal=getattr(t,cutVar)[varIdx]                        
                if cutVal<self.cuts[cutKey][0] or cutVal>=self.cuts[cutKey][1]:  passExtraCuts=False
            if not passExtraCuts: self.rec_passSel[iSystVar] = False
                

            if self.rec_passSel[iSystVar]:

                #select particles
                selP4=[]
                for n in xrange(0,t.n):

                    #kinematics cuts
                    if t.pt[n]<self.ptthreshold[1] : continue
                    if abs(t.eta[n])>self.etathreshold : continue

                    #check if particle passes slice cuts
                    failSliceCuts=False
                    for cutKey in self.cuts:
                        if cutKey!='region': continue
                        cutVar,cutIdx=self.cuts[cutKey]
                        phirec=getattr(t,AXISANGLE[cutVar])[varIdx]
                        idxrec=self.getRegionFor( ROOT.TVector2.Phi_mpi_pi(t.phi[n]-phirec) )
                        if idxrec != cutIdx : failSliceCuts=True
                    if failSliceCuts: continue

                    #apply a tracking efficiency scale factors
                    if isMC:
                        tkSF=1.0
                        if isMC:
                            if varyTkEff>0:
                                if varyTkEff==1 : tkSF=self.tkEffSF[mcera]['vtx'].Eval(min(t.nvtx,40))**2
                                if varyTkEff==2 : tkSF=self.tkEffSF['BCDEF']['vtx'].Eval(min(t.nvtx,40))**2
                                if varyTkEff==3 : tkSF=self.tkEffSF['GH']['vtx'].Eval(min(t.nvtx,40))**2
                                if varyTkEff==4 : tkSF=self.tkEffSF[mcera]['abseta'].Eval(abs(t.eta[n]))**2
                                if varyTkEff==5 : tkSF=self.tkEffSF[mcera+'_dstar']['eta'].Eval(t.eta[n])**2

                        #attempt to remove track if scale factor smaller than one
                        if tkSF<1.0:
                            rnd=ROOT.gRandom.Uniform(1.0)                        
                            if rnd>tkSF : continue

                    #save charged particle
                    p4.SetPtEtaPhiM(t.pt[n],t.eta[n],t.phi[n],self.piMass)
                    selP4.append(ROOT.TLorentzVector(p4))

                    self.rec_chmult[iSystVar]  +=1
                    self.rec_chflux[iSystVar]  += p4.Pt()
                    self.rec_chfluxz[iSystVar] += abs(p4.Pz())

                #loop over gen particles and bring unmatched to reco-level if tracking scale factor >1
                if isMC:
                    for n in xrange(0,t.gen_n):

                        #require non-matched
                        if t.gen_rec[n]>=0 : continue

                        #require pass kinematics
                        if t.gen_pt[n] < self.ptthreshold[0] : continue
                        if abs(t.gen_eta[n])>self.etathreshold : continue

                        #check if it passes slice cuts
                        failSliceCuts=False
                        for cutKey in self.cuts:
                            if cutKey!='region':
                                cutVar=cutKey
                                genCutVar='gen_'+cutVar
                                cutVal=getattr(t,genCutVar)
                                if cutVal<self.cuts[cutKey][0] or cutVal>=self.cuts[cutKey][1]:  failSliceCuts=True
                            else:
                                cutVar,cutIdx=self.cuts[cutKey]
                                genCutVar='gen_'+cutVar
                                phigen=getattr(t,'gen_'+AXISANGLE[cutVar])
                                idxgen=self.getRegionFor( ROOT.TVector2.Phi_mpi_pi(t.gen_phi[n]-phigen) )
                                if idxgen != cutIdx : failSliceCuts=True
                        if failSliceCuts: continue

                        #tracking scale factor
                        tkSF=1.0
                        if varyTkEff==1 : tkSF=self.tkEffSF[mcera]['vtx'].Eval(min(t.nvtx,40))**2
                        if varyTkEff==2 : tkSF=self.tkEffSF['BCDEF']['vtx'].Eval(min(t.nvtx,40))**2
                        if varyTkEff==3 : tkSF=self.tkEffSF['GH']['vtx'].Eval(min(t.nvtx,40))**2
                        if varyTkEff==4 : tkSF=self.tkEffSF[mcera]['abseta'].Eval(abs(t.gen_eta[n]))**2
                        if varyTkEff==5 : tkSF=self.tkEffSF[mcera+'_dstar']['eta'].Eval(t.gen_eta[n])**2

                        #promote gen particle as reco candidate
                        if tkSF<=1.0: continue
                        rnd=ROOT.gRandom.Uniform(1.0)  
                        if rnd>(tkSF-1.0) : continue
                        p4.SetPtEtaPhiM(t.gen_pt[n],t.gen_eta[n],t.gen_phi[n],self.piMass)
                        selP4.append(ROOT.TLorentzVector(p4))

                        self.rec_chmult[iSystVar]  +=1
                        self.rec_chflux[iSystVar]  += p4.Pt()
                        self.rec_chfluxz[iSystVar] += abs(p4.Pz())

                #average pt/pz
                self.rec_chavgpt[iSystVar] = self.rec_chflux[iSystVar]/self.rec_chmult[iSystVar]  if self.rec_chmult[iSystVar]>0 else 0.
                self.rec_chavgpz[iSystVar] = self.rec_chfluxz[iSystVar]/self.rec_chmult[iSystVar] if self.rec_chmult[iSystVar]>0 else 0.
                
                #event shapes
                self.evshapes.analyseNewEvent(selP4)
                self.rec_sphericity[iSystVar]     = self.evshapes.sphericity
                self.rec_aplanarity[iSystVar]     = self.evshapes.aplanarity
                self.rec_C[iSystVar]              = self.evshapes.C
                self.rec_D[iSystVar]              = self.evshapes.D
                
        #gen level (only one selection variation needed)
        self.gen_passSel=(t.gen_passSel&0x1)

        passExtraCuts=True
        for cutKey in self.cuts:
            if cutKey=='region': continue
            cutVar=cutKey
            genCutVar='gen_'+cutVar
            cutVal=getattr(t,genCutVar)
            if cutVal<self.cuts[cutKey][0] or cutVal>=self.cuts[cutKey][1]:  passExtraCuts=False
        if not passExtraCuts: self.gen_passSel=False

        if self.gen_passSel:

            selP4=[]
            for n in xrange(0,t.gen_n):

                #kinematics cuts
                if t.gen_pt[n] < self.ptthreshold[0] : continue
                if abs(t.gen_eta[n])>self.etathreshold : continue

                #check if it passes slice cuts
                failSliceCuts=False
                for cutKey in self.cuts:
                    if cutKey!='region': continue
                    cutVar,cutIdx=self.cuts[cutKey]
                    genCutVar='gen_'+cutVar
                    phigen=getattr(t,'gen_'+AXISANGLE[cutVar])
                    idxgen=self.getRegionFor( ROOT.TVector2.Phi_mpi_pi(t.gen_phi[n]-phigen) )
                    if idxgen != cutIdx : failSliceCuts=True
                if failSliceCuts: continue

                p4.SetPtEtaPhiM(t.gen_pt[n],t.gen_eta[n],t.gen_phi[n],self.piMass)
                selP4.append(ROOT.TLorentzVector(p4))

                self.gen_chmult +=1
                self.gen_chflux += p4.Pt()
                self.gen_chfluxz += abs(p4.Pz())

            #average pt/pz
            self.gen_chavgpt = self.gen_chflux/self.gen_chmult  if self.gen_chmult>0 else 0.
            self.gen_chavgpz = self.gen_chfluxz/self.gen_chmult if self.gen_chmult>0 else 0.

            #event shapes
            self.evshapes.analyseNewEvent(selP4)
            self.gen_sphericity     = self.evshapes.sphericity
            self.gen_aplanarity     = self.evshapes.aplanarity
            self.gen_C              = self.evshapes.C
            self.gen_D              = self.evshapes.D

        #print event if required
        if debug : self.show()

