#ifndef _BjetChargeTreeProducer_h_
#define _BjetChargeTreeProducer_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

struct BJetSummary_t
{
  float pt,eta,phi,m,csv;
  int nch,ch[200];
  float chpt[200], cheta[200], chphi[200], chm[200]; 
  float ptD,ptDs,width,tau21,tau32,tau43,zg;
  float ch_05,ch_08,ch_1,ch_15;
  float ch2_05,ch2_08,ch2_1,ch2_15;
  float ch3_05,ch3_08,ch3_1,ch3_15;

  float g_pt,g_eta,g_phi,g_m;
  int g_bId,g_pId;
  float g_xb;
};


void RunBjetChargeTreeProducer(TString filename,
                               TString outname,
                               Bool_t debug=false);

void createBJetSummaryTree(TTree *t,BJetSummary_t &summary);


#endif
