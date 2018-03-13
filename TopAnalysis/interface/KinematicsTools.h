#ifndef _kinematicstools_h_
#define _kinematicstools_h_

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

struct JetPullInfo_t
{
  Int_t n,nch;
  TVector2 pull,chPull;
};
JetPullInfo_t getPullVector( MiniEvent_t &ev, int ijet);

Float_t computeMT(TLorentzVector &a, TLorentzVector &b);

#endif
