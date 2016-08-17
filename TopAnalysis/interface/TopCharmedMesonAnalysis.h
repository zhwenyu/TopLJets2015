#ifndef _TopCharmedMesonAnalysis_h_
#define _TopCharmedMesonAnalysis_h_

#include "TString.h"
#include "TH1F.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"

void RunTopCharmedMesonAnalysis(TString filename,
				TString outname,
				Int_t channelSelection, 
				Int_t chargeSelection, 
				FlavourSplitting flavourSplitting,
				TH1F *normH, 
				Bool_t runSysts,
				TString era,
				Bool_t debug=false);


/**
   @short summarizes the information on a jet needed for the charmed meson analysis
 */
typedef std::pair<TLorentzVector,int> IdTrack;
class Jet {

 public:
  
  Jet(TLorentzVector p4, float csv, int idx);
  ~Jet();
  void addTrack(TLorentzVector p4, int pfid);
  TLorentzVector &getVec();
  float &getCSV();
  int &getJetIndex();
  std::vector<IdTrack> &getTracks();
  IdTrack getTrack(int idx,float mass);
  void sortTracksByPt();
 
 private:
  TLorentzVector p4_;
  float csv_;
  int idx_;
  std::vector<IdTrack> trks_;
};

bool sortJetsByPt(Jet i, Jet j);
bool sortJetsByCSV(Jet i, Jet j);
bool sortIdTracksByPt(IdTrack i, IdTrack j);

#endif
