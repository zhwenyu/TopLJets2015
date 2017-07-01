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
  float ptD, ptDs, width, tau21, tau32, tau43, zg;

  float g_pt,g_eta,g_phi,g_m;
  int g_bHad,g_pId;
  float g_xb;
  int g_nch,g_ch[200];
  float g_chpt[200], g_cheta[200], g_chphi[200], g_chm[200]; 
  float g_ptD, g_ptDs, g_width, g_tau21, g_tau32, g_tau43, g_zg;
};


void RunBjetChargeTreeProducer(TString filename,
                               TString outname,
                               Bool_t debug=false);

void createBJetSummaryTree(TTree *t,BJetSummary_t &summary);


#endif
