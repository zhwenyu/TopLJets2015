#ifndef _BjetChargeTreeProducer_h_
#define _BjetChargeTreeProducer_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

struct BJetSummary_t
{
  float pt,eta,phi,m;
  float vtxmass,vtxchi2,vtxL3d,vtxL3dSig;
  int nch,vtxnch,nmu,vtxnmu;
  float ch,vtxch,much,vtxmuch;

  float g_pt,g_eta,g_phi,g_m;
  int g_bId,g_pId;
  float g_xb;
};


void RunBjetChargeTreeProducer(TString filename,
                               TString outname,
                               Bool_t debug=false);

void createBJetSummaryTree(TTree *t,BJetSummary_t &summary);


#endif
