#ifndef _correction_tools_h_
#define _correction_tools_h_

#include <vector>
#include <string>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>

#include "TGraphAsymmErrors.h"

typedef std::pair<TString,float> RunPeriod_t;
std::vector<RunPeriod_t> getRunPeriods(TString era);
TString assignRunPeriod(std::vector<RunPeriod_t> &runPeriods,TRandom *rand=0);

std::pair<std::map<Int_t,Float_t>, TH1F *> parseLumiInfo(TString era);
std::map<Int_t,Float_t> lumiPerRun(TString era="era2017");

//pileup weighting
std::vector<TGraph *> getPileupWeights(TString era, TH1 *genPU, TString period = "");
std::map<TString, std::vector<TGraph *> > getPileupWeightsMap(TString era, TH1 *genPU);

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

//apply tracking efficiency
std::map<TString, std::map<TString, std::vector<double> > > getTrackingEfficiencyMap(TString era);
void applyTrackingEfficiencySF(MiniEvent_t &ev, double sf, double minEta = -2.4, double maxEta = 2.4);
void applyEtaDepTrackingEfficiencySF(MiniEvent_t &ev, std::vector<double> sfs, std::vector<double> etaBins);

#endif
