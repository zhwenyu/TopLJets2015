#ifndef _jec_tools_h_
#define _jec_tools_h_

#include <vector>
#include <string>
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "TGraphAsymmErrors.h"

class JECTools
{
 public:
  JECTools();

  FactorizedJetCorrector *getFactorizedJetEnergyCorrector(TString,bool);
  std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt);
  //apply jec uncertainty
  void applyJetCorrectionUncertainty(MiniEvent_t &ev, JetCorrectionUncertainty *jecUnc, TString jecVar, TString direction);
  void applyJetCorrectionUncertainty(TLorentzVector &jp4,JetCorrectionUncertainty *jecUnc,TString direction);
  //apply jet energy resolutions
  void smearJetEnergies(MiniEvent_t &ev, std::string option = "central");
  void smearJetEnergies(MiniEvent_t &ev, JME::JetResolution* jer, std::string option = "central");
  void smearJetEnergy(TLorentzVector &jp4, float genJet_pt,std::string option = "central");
  void smearJetEnergyStochastic(TLorentzVector &jp4, TRandom* random, double resolution, std::string option = "central");
  //b fragmentation, see https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer
  std::map<TString, TGraph*> getBFragmentationWeights(TString era);
  double computeBFragmentationWeight(MiniEvent_t &ev, TGraph* wgtGr);
  std::map<TString, std::map<int, double> > getSemilepBRWeights(TString era);
  double computeSemilepBRWeight(MiniEvent_t &ev, std::map<int, double> corr, int pid = 0, bool useabs = true);

};

#endif
