#ifndef _WbChargeAsymmetry_h_
#define _WbChargeAsymmetry_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

void RunWbChargeAsymmetry(TString filename,
                          TString outname,
                          Int_t channelSelection, 
                          Int_t chargeSelection, 
                          TH1F *normH, 
                          TString era,
                          Bool_t debug=false);

struct WbChargeAsymmetryEvent_t
{

  UInt_t run,event,lumi,cat, nvtx;
  Int_t nw;
  Float_t weight[1000];

  Bool_t reco_sel, gen_sel;

  Int_t l_id,l_c;  
  Float_t l_pt, l_eta, l_phi, l_m;
  Float_t j_pt, j_eta, j_phi, j_m, j_csv;
  Float_t j_vtxmass,j_vtx3DVal,j_vtx3DSig,j_vtxpx,j_vtxpy,j_vtxpz;
  Int_t j_vtxntk;
  Int_t ntk,tk_c[100], tk_id[100];
  Float_t tk_pt[100],tk_eta[100],tk_phi[100];
  Int_t nds,ds_spc[100];
  Float_t ds_m12[100],ds_dm[100],ds_pt[100],ds_eta[100],ds_phi[100];
  Float_t met_pt,met_phi;

  Int_t nj,ngj;

  Int_t gl_id,gl_c;
  Float_t gl_pt, gl_eta, gl_phi, gl_m;  
  Int_t gj_flavor,gj_bid;
  Bool_t gj_isSemiLepBhad;
  Float_t gj_pt, gj_eta, gj_phi, gj_m,gj_xb;
  Int_t ngtk,gtk_c[100],gtk_id[100];
  Float_t gtk_pt[100],gtk_eta[100],gtk_phi[100];

};

struct DsCand_t
{
  TLorentzVector pi,kaon,softpi;
  int softpicharge;
  float mass12,dM;
};


void createWbChargeAsymmetryEventTree(TTree *t,WbChargeAsymmetryEvent_t &tjsev);
void resetWbChargeAsymmetryEvent(WbChargeAsymmetryEvent_t &tjsev);
std::vector<DsCand_t> searchForDsTag(std::vector<Particle> &pColl);

#endif
