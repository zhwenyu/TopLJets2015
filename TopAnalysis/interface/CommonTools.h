#ifndef _common_tools_h_
#define _common_tools_h_

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

#include "TVector2.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TRandom.h"

#include <vector>

Float_t computeMT(TLorentzVector &a, TLorentzVector &b);
FactorizedJetCorrector *getFactorizedJetEnergyCorrector(TString,bool);
std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt);
float getLeptonEnergyScaleUncertainty(int l_id,float l_pt,float l_eta);

struct JetPullInfo_t
{
  Int_t n,nch;
  TVector2 pull,chPull;
};
JetPullInfo_t getPullVector( MiniEvent_t &ev, int ijet);

class HistTool {

 public:
  HistTool(int nsyst = 20);
  ~HistTool() {}
  
  void setNsyst(int nsyst) { nsyst_ = nsyst; }
  void addHist(TString title, TH1* hist);
  void fill(TString title, double value, std::vector<double> weights);
  std::map<TString, TH1 *> &getPlots()   { return allPlots_; }
  std::map<TString, TH2 *> &get2dPlots() { return all2dPlots_; }
  
 private:
  int nsyst_;
  std::map<TString, TH1 *> allPlots_;
  std::map<TString, TH2 *> all2dPlots_;

};

#endif
