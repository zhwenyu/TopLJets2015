#ifndef _tophiforest_h_
#define _tophiforest_h_

#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"

struct LJEvent_t
{
  Float_t w;
  ULong64_t run,lumi,event;
  Float_t hiHFplus,hiHFminus,hiHFplusEta4,hiHFminusEta4;
  Int_t nj,ngj,ngp,nb,npt;
  Int_t j_btag[20];  
  Int_t pt_id[20],pt_pdgid[20];
  Float_t pt_pt[20],pt_eta[20],pt_phi[20],pt_m[20];
  Int_t gp_pdgid[20];
  Float_t gp_pt[20],gp_eta[20],gp_phi[20],gp_m[20];
  Float_t j_pt[20],j_eta[20],j_phi[20],j_m[20], j_area[20];
  Float_t j_PfCHF[20],j_PfNHF[20],j_PfCEF[20],j_PfNEF[20], j_PfMUF[20],j_PfCHM[20];//jet id
  Int_t j_PfNHM[20], j_PfCEM[20],j_PfNEM[20], j_PfMUM[20]; //jet id
  Float_t j_refpt[20],j_refarea[20],j_refdr[20], j_refparton_flavorForB[20], j_refparton_flavor[20];
  Float_t gj_pt[20],gj_eta[20],gj_phi[20],gj_m[20],gj_dr[20],gj_index[20];
  Int_t l_id;
  Float_t l_pt, l_eta, l_phi, l_m;
  Float_t gl_pt, gl_eta, gl_phi, gl_m;
  Float_t met_pt, met_phi;
  Int_t ntracks,ntracks_hp;
};

bool sortDijetKeys(std::pair< std::pair<int,int>, float > a,
		   std::pair< std::pair<int,int>, float > b);
std::pair<int,int> getDijetsSystemCandidate(std::vector<TLorentzVector> &lightJets);
void defineTreeBranches(TTree *,LJEvent_t &,bool);
void RunTop16023(TString filename,
		     TString outname,
		     Int_t channelSelection, 
		     Int_t chargeSelection, 
		     TH1F *normH, 
		     Bool_t runSysts,
		     TString era);
void RunHin17002(TString filename,
                 TString outname,
                 Int_t channelSelection, 
                 Int_t chargeSelection, 
                 TH1F *normH, 
                 Bool_t runSysts,
                 TString era);

#endif
