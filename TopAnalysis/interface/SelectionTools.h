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

  enum AnalysisType	{TOP=0,VBF=1};

  SelectionTool(TString dataset="",bool debug=false,TH1 *triggerList=0, AnalysisType analysisType = TOP);
  ~SelectionTool() {}

  enum FlavourSplitting {NOFLAVOURSPLITTING=0, UDSGSPLITTING=1, CSPLITTING=4, BSPLITTING=5 };
  enum QualityFlags     {VETO, LOOSE, MEDIUM, TIGHT, CONTROL, FAKE};



  //
  //RECO LEVEL SELECTORS
  // 
  TString flagFinalState(MiniEvent_t &ev, std::vector<Particle> leptons={}, std::vector<Particle> photons={}, bool isCR=false); 
  std::vector<Particle> leptons_,vetoLeptons_,photons_;
  std::vector<Jet> jets_;
  TLorentzVector met_;
  bool hasTriggerBit(TString triggerName,unsigned int word);
  std::vector<Particle> flaggedLeptons(MiniEvent_t &ev);
  std::vector<Particle> selLeptons(std::vector<Particle> &flaggedLeptons,int qualBit=LOOSE, double minPt=0., double maxEta=99., std::vector<Particle> veto={});
  std::vector<Particle> &getSelLeptons()  { return leptons_; }
  std::vector<Particle> &getVetoLeptons() { return vetoLeptons_; }
  std::vector<Particle> flaggedPhotons(MiniEvent_t &ev);
  std::vector<Particle> selPhotons(std::vector<Particle> &flaggedPhotons,int qualBit=LOOSE, std::vector<Particle> leptons = {}, double minPt = 30., double maxEta = 2.5, std::vector<Particle> veto = {});
  std::vector<Particle> &getSelPhotons()  { return photons_; }
  std::vector<Jet>      getGoodJets(MiniEvent_t &ev, double minPt = 30., double maxEta = 4.7, std::vector<Particle> leptons = {},std::vector<Particle> photons = {}); // changed for the moment to VBF
  //void                  getGoodJets(MiniEvent_t &ev, double minPt = 30., double maxEta = 4.7, std::vector<Particle> leptons = {},std::vector<Particle> photons = {}); 
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

  void setPhotonSelection(std::vector<TString> photonTrigs={"HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v","HLT_Photon200_v"},
                          SelectionTool::QualityFlags offlinePhoton=TIGHT, double addedFakeSieie = 0.0027){
    photonTriggers_=photonTrigs;
    offlinePhoton_=offlinePhoton;
    setFakePhotonSelection(addedFakeSieie);
    
  }
  void setFakePhotonSelection(double addedsieie = 0.0027){
    /////////////////////////////////
    // Modified Loose Requirements //
    // Only sieie has been changed //
    /////////////////////////////////
    std::vector<double> eeloose, ebloose;

    ebloose.push_back(0.105);
    ebloose.push_back(0.0103+addedsieie);
    eeloose.push_back(0.029);
    eeloose.push_back(0.0276+addedsieie);
    //charged iso
    ebloose.push_back(2.839);
    eeloose.push_back(2.150);
    //neutral iso
    ebloose.push_back(9.188);
    ebloose.push_back(0.0126);
    ebloose.push_back(0.000026);
    eeloose.push_back(10.471);
    eeloose.push_back(0.0119);
    eeloose.push_back(0.000025);
    //photoniso
    ebloose.push_back(2.956);
    ebloose.push_back(0.0035);
    eeloose.push_back(4.895);
    eeloose.push_back(0.0040);
    
    fakeIdCuts["EB"] = ebloose;
    fakeIdCuts["EE"] = eeloose;
  }

  bool isInclusivePhoton(MiniEvent_t ev, int idx){
    // from AN-14-242
    double pt(ev.gamma_pt[idx]);
    double eta(fabs(ev.gamma_eta[idx]));
    TString region = "EB";
    if(fabs(eta) > 1.4442)
      region = "EE";

    double cutNeutIso   = (fakeIdCuts[region][3] + pt*fakeIdCuts[region][4] + pt*pt*fakeIdCuts[region][5]);
    double cutGIso      = (fakeIdCuts[region][6] + pt*fakeIdCuts[region][7]);
    bool HoE     = ev.gamma_hoe[idx]                    < fakeIdCuts[region][0];
    bool sieie   = ev.gamma_sieie[idx]                  < fakeIdCuts[region][1];
    bool chIso   = ev.gamma_chargedHadronIso[idx]       < std::min(5* fakeIdCuts[region][2], 0.2*pt);
    bool neutIso = ev.gamma_neutralHadronIso[idx]       < std::min(5* cutNeutIso           , 0.2*pt);
    bool gIso    = ev.gamma_photonIso[idx]              < std::min(5* cutGIso              , 0.2*pt);

    bool ret = (HoE && sieie && chIso && neutIso && gIso);
    return ret;
    
  }
  bool isFakePhoton(MiniEvent_t ev, int idx){

    // from AN-14-242
    double pt(ev.gamma_pt[idx]);
    double eta(fabs(ev.gamma_eta[idx]));
    TString region = "EB";
    if(fabs(eta) > 1.4442)
      region = "EE";

    double cutNeutIso   = (fakeIdCuts[region][3] + pt*fakeIdCuts[region][4] + pt*pt*fakeIdCuts[region][5]);
    double cutGIso      = (fakeIdCuts[region][6] + pt*fakeIdCuts[region][7]);

    double looseChIso   = ev.gamma_chargedHadronIso[idx] < fakeIdCuts[region][2];
    double looseNeutIso = ev.gamma_neutralHadronIso[idx] < cutNeutIso;
    double looseGIso    = ev.gamma_photonIso[idx]        < cutGIso;

    bool notLooseIso = !(looseChIso || looseNeutIso || looseGIso);

    bool ret = (notLooseIso && this->isInclusivePhoton(ev,idx));
    return ret;
  }

 private:
  bool debug_;
  AnalysisType anType_;
  bool isSingleElectronPD_,isSingleMuonPD_,isDoubleEGPD_,isDoubleMuonPD_,isMuonEGPD_,isPhotonPD_,isJetHTPD_;
  std::map<TString,unsigned int> triggerBits_;
  std::vector<TString> photonTriggers_;
  int offlinePhoton_;
  //{[EB,EE], hoe, sieie, chiso, nuetiso, photonsio}
  std::map<TString, std::vector<double> > fakeIdCuts;
};

#endif
