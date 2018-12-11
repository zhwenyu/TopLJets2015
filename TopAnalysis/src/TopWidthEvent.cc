#include "TopLJets2015/TopAnalysis/interface/TopWidthEvent.h"
#include <algorithm>
#include <tuple>
#include <iostream>

using namespace std;

//
TopWidthEvent::TopWidthEvent(std::vector<Particle> &leptons,std::vector<Jet> &jets) : cat(""), dilcode(0) {

  initSelectionCuts();

  if(leptons.size()<2 || jets.size()<2) return;
  
  //lepton selection
  if(leptons[0].Pt()<leadLeptonPt_    || fabs(leptons[0].Eta())>maxLeptonEta_) return;
  if(leptons[1].Pt()<subLeadLeptonPt_ || fabs(leptons[1].Eta())>maxLeptonEta_) return;
  l1pt  = leptons[0].Pt();
  l1eta = leptons[0].Eta();
  l2pt  = leptons[1].Pt();
  l2eta = leptons[1].Eta();
  
  //dilepton selection
  if(leptons[0].charge()*leptons[1].charge()>0) return;
  mll   = ( (leptons[0]+leptons[1]).M() );
  ptll  = ( (leptons[0]+leptons[1]).Pt() );
  if(mll<20) return;
  bool isZ( fabs(mll-91)<15 );
  
  //jet selection
  std::vector<Jet> bJets,lJets;
  for(auto j : jets) {
    if(j.Pt()<jetPt_) continue;
    if(fabs(j.Eta())>maxJetEta_) continue;    
    bool isBTagged(j.flavor()==5);
    if(isBTagged) bJets.push_back(j);
    else          lJets.push_back(j);
  }
  nbjets=bJets.size();
  njets=nbjets+lJets.size();
  if(njets<2) return;
  if(nbjets<1) return;

  j1pt  = bJets[0].Pt();
  j1eta = bJets[0].Eta();
  j2pt  = bJets.size()>1 ? bJets[1].Pt()  : lJets[0].Pt();
  j2eta = bJets.size()>1 ? bJets[1].Eta() : lJets[0].Eta();
  
  //lepton-b pairing
  if(bJets.size()==1) {
    
    //choose pair with leading ptlb
    TLorentzVector lb1(leptons[0]+bJets[0]),lb2(leptons[1]+bJets[0]);
    lbPairs_.push_back( lb1.Pt()>lb2.Pt() ? lb1 : lb2 );
  }
  else {
    
    //sort pairs by decreasing ptlb
    std::vector< std::tuple<int,int,TLorentzVector,float> > allPairs;
    for(size_t il=0; il<2; il++)
      for(size_t ib=0; ib<2; ib++) {
        TLorentzVector lb(leptons[il]+bJets[il]);
        allPairs.push_back( std::make_tuple(il,ib,lb,lb.Pt()) );
      }
    std::sort(allPairs.begin(),allPairs.end(),
              [](const tuple<int,int,TLorentzVector,float>& a,
                 const tuple<int,int,TLorentzVector,float>& b) -> bool
              {
                return std::get<3>(a) > std::get<3>(b);
              });

    lbPairs_.push_back( std::get<2>(allPairs[0]) );
    for(size_t ipair=1; ipair<4; ipair++) {
      if( std::get<0>(allPairs[0]) ==  std::get<0>(allPairs[ipair]) ) continue;
      if( std::get<0>(allPairs[1]) ==  std::get<1>(allPairs[ipair]) ) continue;
      lbPairs_.push_back( std::get<2>(allPairs[ipair]) );
      break;
    }
  }


  //assign category
  dilcode=(abs(leptons[0].id()*leptons[1].id()));
  if(dilcode==11*11) cat=isZ ? "zee":"ee";
  if(dilcode==13*13) cat=isZ ? "zmm":"mm";
  if(dilcode==11*13) cat="em";
}


//
void TopWidthEvent::initSelectionCuts() {
  leadLeptonPt_    = 30.;
  subLeadLeptonPt_ = 20.;
  maxLeptonEta_    = 2.4;
  zMassWindow_     = 15.;
  minMll_          = 20.;
  jetPt_           = 30.;
  maxJetEta_       = 2.4;
} 
