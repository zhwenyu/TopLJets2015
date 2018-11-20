#ifndef _VBFDISCRIMINATORINPUTS_H_
#define _VBFDISCRIMINATORINPUTS_H_

#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

namespace vbf{

  class DiscriminatorInputs{
  public:
    float leadj_pt,leadj_gawidth,leadj_c2_02,lead_qg;
    float subleadj_pt,subleadj_gawidth,subleadj_c2_02,sublead_qg;
    float mjj,detajj,dphijj,jjpt,jjetas,ystar,balance,relbpt,dphibjj,dphivj0,dphivj1;
    float isotropy,circularity,sphericity,aplanarity,C,D;
    float scalarht,mht;
  DiscriminatorInputs(): 
    leadj_pt(0),leadj_gawidth(-99),leadj_c2_02(-99),lead_qg(-99),
      subleadj_pt(0),subleadj_gawidth(-99),subleadj_c2_02(-99),sublead_qg(-99),
      mjj(0),detajj(-99),dphijj(-99),jjpt(0),jjetas(-99),ystar(-99),balance(-99),relbpt(-99),dphibjj(-99),dphivj0(-99),dphivj1(-99),
      isotropy(-99),circularity(-99),sphericity(-99),aplanarity(-99),C(-99),D(-99),scalarht(0),mht(0) {
    }
  DiscriminatorInputs(const DiscriminatorInputs &o): 
    leadj_pt(o.leadj_pt),leadj_gawidth(o.leadj_gawidth),leadj_c2_02(o.leadj_c2_02),lead_qg(o.leadj_c2_02),
      subleadj_pt(o.subleadj_pt),subleadj_gawidth(o.subleadj_pt),subleadj_c2_02(o.subleadj_pt),sublead_qg(o.subleadj_pt),
      mjj(o.mjj),detajj(o.detajj),dphijj(o.dphijj),jjpt(o.jjpt),jjetas(o.jjetas),ystar(o.ystar),balance(o.balance),relbpt(o.relbpt),dphibjj(o.dphibjj),dphivj0(o.dphivj0),dphivj1(o.dphivj0),
      isotropy(o.isotropy),circularity(o.circularity),sphericity(o.sphericity),aplanarity(o.aplanarity),C(o.C),D(o.D),
      scalarht(o.scalarht), mht(o.mht){
    }
    DiscriminatorInputs& operator=(const DiscriminatorInputs& o)
      {
        if(&o == this) return *this;
        leadj_pt=o.leadj_pt;  leadj_gawidth=o.leadj_gawidth;  leadj_c2_02=o.leadj_c2_02;  lead_qg=o.leadj_c2_02;  
        subleadj_pt=o.subleadj_pt;  subleadj_gawidth=o.subleadj_pt;  subleadj_c2_02=o.subleadj_pt;  sublead_qg=o.subleadj_pt;  
        mjj=o.mjj;  detajj=o.detajj;  dphijj=o.dphijj;  jjpt=o.jjpt;  jjetas=o.jjetas;  ystar=o.ystar;  balance=o.balance;  relbpt=o.relbpt;  dphibjj=o.dphibjj;  dphivj=o.dphivj0;  dphivj1=o.dphivj0;  
        isotropy=o.isotropy;  circularity=o.circularity;  sphericity=o.sphericity;  aplanarity=o.aplanarity;  C=o.C;  D=o.D;
        scalarht=o.scalarht; mht=o.mht;
        return *this;
      }
    inline void fillDiscriminatorVariables(Particle &boson,std::vector<Jet> &jets,MiniEvent_t &ev) {

      if(jets.size()>0){
        vars.leadj_pt      = jets[0].Pt();
        vars.leadj_gawidth = ev.j_gawidth[jets[0].getJetIndex()];
        vars.leadj_c2_02   = ev.j_c2_02[jets[0].getJetIndex()];
        vars.lead_qg       = ev.j_qg[jets[0].getJetIndex()];
        vars.dphivj0       = fabs(jets[0].DeltaPhi(boson));
        vars.centraleta    = jets[0].Eta();
      }
    
      if(jets.size()>=2){
      
        vars.subleadj_pt      = jets[1].Pt();
        vars.subleadj_gawidth = ev.j_gawidth[jets[1].getJetIndex()];
        vars.subleadj_c2_02   = ev.j_c2_02[jets[1].getJetIndex()];
        vars.subleadj_gq      = ev.j_qg[jets[0].getJetIndex()];
        vars.dphivj1          = fabs(jets[1].DeltaPhi(boson));
        vars.centraleta       = min(fabs(jets[0].Eta()),fabs(jets[1].Eta()));
        vars.forwardeta       = max(fabs(jets[0].Eta()),fabs(jets[1].Eta()));

        TLorentzVector jj(jets[0]+jets[1]);
        vars.mjj     = jj.M();
        vars.detajj  = fabs(jets[0].Eta()-jets[1].Eta());
        vars.dphijj  = jets[0].DeltaPhi(jets[1]);
        vars.jjpt    = (jets[0]+jets[1]).Pt();
        vars.jjetas  = jets[0].Eta()*jets[1].Eta();
        vars.ystar   = boson.Rapidity()-0.5*(jets[0].Rapidity()+jets[1].Rapidity());
        vars.balance = (boson+jets[0]+jets[1]).Pt();
        vars.relbpt  = (jets[0].Pt()+jets[1].Pt())/boson.Pt();
        vars.dphibjj = boson.DeltaPhi( jets[0]+jets[1] );    
      }


      scalarht = 0.;
      TLorentzVector mhtP4(0,0,0,0);
      for(auto j : jets) {
        scalarht += j.Pt();
        mhtP4 += j;
      }
      mht = mhtP4.Pt();

    
      //event shapes
      std::vector<math::XYZVector> inputVectors;
      inputVectors.push_back( math::XYZVector(boson.Px(),boson.Py(),boson.Pz()) );
      for(size_t ij=0; ij<min(size_t(2),jets.size());ij++) {
        inputVectors.push_back( math::XYZVector(jets[ij].Px(),jets[ij].Py(),jets[ij].Pz()) );
      }
      EventShapeVariables esv(inputVectors);
      vars.isotropy    = esv.isotropy();
      vars.circularity = esv.circularity();
      vars.sphericity  = esv.sphericity(1.);
      vars.aplanarity  = esv.aplanarity(1.);
      vars.C           = esv.C(1.);
      vars.D           = esv.D(1.);
    }
  };

}

#endif
