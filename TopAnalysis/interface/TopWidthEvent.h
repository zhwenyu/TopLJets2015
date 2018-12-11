#ifndef _topwidthevent_h_
#define _topwidthevent_h_

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"

class TopWidthEvent {

 public:

  TopWidthEvent(std::vector<Particle> &leptons,std::vector<Jet> &jets);
  std::vector<TLorentzVector> getPairs() { return lbPairs_; }
  ~TopWidthEvent() { }

  //access directly
  TString cat;
  unsigned int dilcode;
  float mll,ptll,l1pt,l1eta,l2pt,l2eta,njets,nbjets,j1pt,j1eta,j2pt,j2eta;

 private:
  std::vector<TLorentzVector> lbPairs_;
  void initSelectionCuts();
  float leadLeptonPt_,subLeadLeptonPt_,maxLeptonEta_,zMassWindow_,minMll_,jetPt_,maxJetEta_;

};


#endif
