#ifndef _LeptonEfficiencyWrapper_h_
#define _LeptonEfficiencyWrapper_h_

#include "TString.h"
#include "TH2F.h"
#include "TLorentzVector.h"

#include <vector>
#include <map>

typedef std::map<TString,TH2 *> LeptonEfficiency_t;
typedef std::pair<float,float> EffCorrection_t;

class LeptonEfficiencyWrapper : public LeptonEfficiency_t
{
 public:
  LeptonEfficiencyWrapper(bool isData,TString era);
  EffCorrection_t getTriggerCorrection(std::vector<int> pdgId,std::vector<TLorentzVector> leptons);
  EffCorrection_t getOfflineCorrection(int pdgId,float pt,float eta);
  ~LeptonEfficiencyWrapper();  
 private:
  int era_;
  void init(TString era);
};

#endif
