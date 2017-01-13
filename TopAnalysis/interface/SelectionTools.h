#ifndef _selection_tools_h_
#define _selection_tools_h_

#include <vector>

#include "TString.h"

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"

class SelectionTool {

 public:
  SelectionTool(TString dataset="",bool debug=false);
  ~SelectionTool() {}

  enum FlavourSplitting {NOFLAVOURSPLITTING=0, UDSGSPLITTING=1, CSPLITTING=4, BSPLITTING=5 };
  enum TOPLeptonQualityFlags { PASSLLID=0, PASSLID, PASSLVETO, PASSLIDNONISO }; 

  //this function will make the standard top selection
  TString flagFinalState(MiniEvent_t &ev, std::vector<Particle> preselleptons={});
  std::vector<Particle> &getSelLeptons()  { return leptons_; }
  std::vector<Particle> &getVetoLeptons() { return vetoLeptons_; }
  std::vector<Jet>      &getJets()        { return jets_; }

  //selection at particle level
  TString flagGenFinalState(MiniEvent_t &ev, std::vector<Particle> preselleptons={});
  std::vector<Particle> &getGenLeptons()  { return genLeptons_; }
  std::vector<Jet>      &getGenJets()     { return genJets_; }
 
 //object selection can also be customized with the following functions
  std::vector<Particle> getTopFlaggedLeptons(MiniEvent_t &ev);
  std::vector<Particle> getGoodLeptons(std::vector<Particle> &topFlaggedLeptons,int qualBit=PASSLLID, double minPt=0., double maxEta=99., std::vector<Particle> *vetoParticles=0);
  std::vector<Jet> getGoodJets(MiniEvent_t &ev, double minPt = 30., double maxEta = 2.4, std::vector<Particle> leptons = {});

  //gen level selection customization
  std::vector<Particle> getGenLeptons(MiniEvent_t &ev, double minPt = 30., double maxEta = 2.1);
  std::vector<Jet> getGenJets(MiniEvent_t &ev, double minPt = 30., double maxEta = 2.4, std::vector<Particle> *leptons = 0);

  //check met filter flags cf. https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
  bool passMETFilters(MiniEvent_t &ev);

 private:
  bool acceptE_, acceptM_, acceptEE_, acceptEM_, acceptMM_;
  std::vector<Particle> leptons_,vetoLeptons_,genLeptons_;
  std::vector<Jet> jets_,genJets_;
  bool debug_;

};

#endif
