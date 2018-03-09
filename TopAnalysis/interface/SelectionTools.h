#ifndef _selection_tools_h_
#define _selection_tools_h_

#include <vector>
#include <map>

#include "TString.h"
#include "TH1.h"

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"

class SelectionTool {

 public:
  SelectionTool(TString dataset="",bool debug=false,TH1 *triggerList=0);
  ~SelectionTool() {}

  enum FlavourSplitting {NOFLAVOURSPLITTING=0, UDSGSPLITTING=1, CSPLITTING=4, BSPLITTING=5 };
  enum QualityFlags     {VETO, LOOSE, MEDIUM, TIGHT, CONTROL};

  //
  //RECO LEVEL SELECTORS
  // 
  TString flagFinalState(MiniEvent_t &ev, std::vector<Particle> leptons={}, std::vector<Particle> photons={});
  std::vector<Particle> leptons_,vetoLeptons_,photons_;
  std::vector<Jet> jets_;
  TLorentzVector met_;
  bool hasTriggerBit(TString triggerName,unsigned int word);
  std::vector<Particle> flaggedLeptons(MiniEvent_t &ev);
  std::vector<Particle> selLeptons(std::vector<Particle> &flaggedLeptons,int qualBit=LOOSE, double minPt=0., double maxEta=99., std::vector<Particle> veto={});
  std::vector<Particle> &getSelLeptons()  { return leptons_; }
  std::vector<Particle> &getVetoLeptons() { return vetoLeptons_; }
  std::vector<Particle> flaggedPhotons(MiniEvent_t &ev);
  std::vector<Particle> selPhotons(std::vector<Particle> &flaggedPhotons,int qualBit=LOOSE, double minPt = 30., double maxEta = 2.5, std::vector<Particle> veto = {});
  std::vector<Particle> &getSelPhotons()  { return photons_; }
  std::vector<Jet>      getGoodJets(MiniEvent_t &ev, double minPt = 30., double maxEta = 2.4, std::vector<Particle> leptons = {},std::vector<Particle> photons = {});
  std::vector<Jet>      &getJets()        { return jets_; }
  bool                   passMETFilters(MiniEvent_t &ev);
  TLorentzVector        &getMET()         { return met_; }

  //
  //PARTICLE LEVEL SELECTORS
  //
  //if no pre-selection has been passed it will trigger the standard pre-selection with the methods below
  TString flagGenFinalState(MiniEvent_t &ev, std::vector<Particle> preselleptons={}, std::vector<Particle> preselphotons={});
  std::vector<Particle> genLeptons_,genPhotons_;
  std::vector<Jet> genJets_;
  std::vector<Particle> getGenLeptons(MiniEvent_t &ev, double minPt = 20., double maxEta = 2.5);
  std::vector<Particle> &getGenLeptons()  { return genLeptons_; }
  std::vector<Particle> getGenPhotons(MiniEvent_t &ev, double minPt = 50., double maxEta = 1.442);
  std::vector<Particle> &getGenPhotons()  { return genPhotons_; }
  std::vector<Jet>      getGenJets(MiniEvent_t &ev, double minPt = 30., double maxEta = 2.4, std::vector<Particle> leptons = {}, std::vector<Particle> photons = {});
  std::vector<Jet>      &getGenJets()     { return genJets_; }

  void setDebug(bool flag) { debug_=flag; }

 private:
  bool debug_;
  bool isSingleElectronPD_,isSingleMuonPD_,isDoubleEGPD_,isDoubleMuonPD_,isMuonEGPD_,isPhotonPD_;
  std::map<TString,unsigned int> triggerBits_;
};

#endif
