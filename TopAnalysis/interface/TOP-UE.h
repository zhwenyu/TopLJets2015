#ifndef _topue_h_
#define _topue_h_

#include "TH1.h"
#include "TString.h"
#include "TFile.h"

#include <map>

struct TopUE_t
{
  Int_t run,event,lumi,cat;
  Int_t nvtx;
  Int_t nw,n,gen_n;
  Int_t passSel,gen_passSel;
  Int_t gen_nj,gen_nb;
  Int_t nj[9],nb[9];
  Int_t isInBFlags[5000];
  Float_t weight[400];
  Float_t gen_pt_ttbar,       gen_phi_ttbar;
  Float_t gen_pt_pseudottbar, gen_phi_pseudottbar;
  Float_t rec_pt_ttbar[9],    rec_phi_ttbar[9];
  Float_t pt[5000],     eta[5000],     phi[5000];
  Float_t gen_pt[5000], gen_eta[5000], gen_phi[5000]; 
  Float_t mll,gen_mll,dphill,gen_dphill;
};

void createTopUETree(TTree *t,TopUE_t &tue);
void resetTopUE(TopUE_t &tue);
void RunTopUE(TString filename,
	      TString outname,
	      Int_t channelSelection, 
	      Int_t chargeSelection, 
	      FlavourSplitting flavourSplitting,
	      TH1F *normH, 
	      Bool_t runSysts,
	      TString era);
#endif
