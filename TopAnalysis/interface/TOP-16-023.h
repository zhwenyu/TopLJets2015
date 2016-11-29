#ifndef _top16023_h_
#define _top16023_h_

#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"

struct LJEvent_t
{
  Float_t w;
  Int_t nj,nb;
  Bool_t j_btag[20];
  Float_t j_pt[20],j_eta[20],j_phi[20],j_m[20];
  Int_t l_id;
  Float_t l_pt, l_eta, l_phi, l_m;
};

void RunTop16023(TString filename,
		     TString outname,
		     Int_t channelSelection, 
		     Int_t chargeSelection, 
		     FlavourSplitting flavourSplitting,
		     TH1F *normH, 
		     Bool_t runSysts,
		     TString era);

#endif
