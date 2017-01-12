#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

SelectionTool::SelectionTool(TString dataset) :
  acceptE_(true), acceptM_(true), acceptEE_(true), acceptEM_(true), acceptMM_(true)
{
  if(dataset.Contains("SingleElectron")) { acceptE_=true;  acceptM_=false; acceptEE_=false;  acceptMM_=false; acceptEM_=false;}
  if(dataset.Contains("SingleMuon"))     { acceptE_=false; acceptM_=true;  acceptEE_=false;  acceptMM_=false; acceptEM_=false;}
  if(dataset.Contains("DoubleEG"))       { acceptE_=false; acceptM_=false; acceptEE_=true;   acceptMM_=false; acceptEM_=false;}
  if(dataset.Contains("DoubleMuon"))     { acceptE_=false; acceptM_=false; acceptEE_=false;  acceptMM_=true;  acceptEM_=false;}
  if(dataset.Contains("MuonEG"))         { acceptE_=false; acceptM_=false; acceptEE_=false;  acceptMM_=false; acceptEM_=true; }
}

//
std::vector<Particle> SelectionTool::getTopFlaggedLeptons(MiniEvent_t ev){
  std::vector<Particle> leptons;
  for (int il=0; il<ev.nl; il++) {

    float pt(ev.l_pt[il]);
    float eta(fabs(ev.l_eta[il]));
    int pid(ev.l_pid[il]);
    float relIso(ev.l_relIso[il]);

    int topLeptonQualityFlagsWord(0);
    if(abs(ev.l_id[il])==11)
      {
	if( pt>20 && eta<2.4 && ((pid>>4) &0x1))                                     topLeptonQualityFlagsWord |= (0x1 << PASSLLID);
	if( pt>30 && eta<2.4 && ((pid>>4) &0x1))                                     topLeptonQualityFlagsWord |= (0x1 << PASSLID);
	if( pt>10 && eta<2.4 && ((pid>>2) &0x1))                                     topLeptonQualityFlagsWord |= (0x1 << PASSLVETO);
	if( pt>26 && eta<2.1 && ((pid>>5) &0x1) && ((pid>>4) &0x1)==0 && relIso>0.4) topLeptonQualityFlagsWord |= (0x1 << PASSLIDNONISO);
      }
    else
      {
	if( pt>20 && eta<2.4 && ((pid>>4) &0x1) && relIso<0.15)  topLeptonQualityFlagsWord |= (0x1 << PASSLLID);
	if( pt>26 && eta<2.1 && ((pid>>4) &0x1) && relIso<0.15)  topLeptonQualityFlagsWord |= (0x1 << PASSLID);
	if( pt>10 && eta<2.4 && ((pid>>1) &0x1) && relIso<0.25)  topLeptonQualityFlagsWord |= (0x1 << PASSLVETO);
	if( pt>26 && eta<2.1 && ((pid>>4) &0x1) && relIso>0.25)  topLeptonQualityFlagsWord |= (0x1 << PASSLIDNONISO);
      }
    if(topLeptonQualityFlagsWord==0) continue;

    TLorentzVector lp4;
    lp4.SetPtEtaPhiM(ev.l_pt[il],ev.l_eta[il],ev.l_phi[il],ev.l_mass[il]);
    leptons.push_back(Particle(lp4, ev.l_charge[il], ev.l_id[il], topLeptonQualityFlagsWord, il, 1.0) );
  }
  
  return leptons;
}

//
std::vector<Particle> SelectionTool::getGoodLeptons(std::vector<Particle> &leptons,int qualBit,double minPt, double maxEta,std::vector<Particle> *vetoParticles){
  std::vector<Particle> selLeptons;

  for(size_t i =0; i<leptons.size(); i++)
    {
      //check quality flag
      if( !leptons[i].hasQualityFlag(qualBit) ) continue;

      //check kinematics
      if(leptons[i].pt()<minPt || fabs(leptons[i].eta())<maxEta) continue;

      //check if this lepton should be vetoed by request
      if(vetoParticles){
	bool skipThisLepton(false);
	for(auto &vetoL : *vetoParticles){
	  if(vetoL.originalReference()!=leptons[i].originalReference()) continue;
	  skipThisLepton=true;
	  break;
	}
	if(skipThisLepton) continue;
      }

      //lepton is selected
      selLeptons.push_back(leptons[i]);
    }

  //all done here
  return selLeptons;
}

//
std::vector<Particle> SelectionTool::getGenLeptons(MiniEvent_t ev, double tightMinPt, double tightMaxEta, bool vetoLoose, double looseMinPt, double looseMaxEta) {
  std::vector<Particle> leptons;
  
  //loop over leptons from pseudotop producer
  for (int i = 0; i < ev.ng; i++) {
    if (abs(ev.g_id[i])>10) {
      bool passTightKin(ev.g_pt[i]>tightMinPt && fabs(ev.g_eta[i])<tightMaxEta);
      bool passLooseKin(ev.g_pt[i]>looseMinPt && fabs(ev.g_eta[i])<looseMaxEta);
      
      if(passTightKin) {
        TLorentzVector lp4;
        lp4.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);
        leptons.push_back( Particle(lp4, -ev.g_id[i]/abs(ev.g_id[i]), ev.g_id[i], 0, 1) );
      }
      else if(vetoLoose && passLooseKin) {
        return {};
      }
    }
  }
  
  return leptons;
}

//
std::vector<Jet> SelectionTool::getGoodJets(MiniEvent_t ev, double minPt, double maxEta, std::vector<Particle> leptons) {
  std::vector<Jet> jets;
  
  for (int k=0; k<ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    //cross clean with leptons
    bool overlapsWithLepton(false);
    for (auto& lepton : leptons) {
      if(jp4.DeltaR(lepton.p4())<0.4) overlapsWithLepton=true;
    }
    if(overlapsWithLepton) continue;
    
    //jet kinematic selection
    if(jp4.Pt() < minPt || abs(jp4.Eta()) > maxEta) continue;

    //flavor based on b tagging
    int flavor = 0;
    if (ev.j_btag[k]) {
      flavor = 5;
    }
    
    Jet jet(jp4, flavor, k);
    
    //fill jet constituents
    for (int p = 0; p < ev.npf; p++) {
      if (ev.pf_j[p] == k) {
        TLorentzVector pp4;
        pp4.SetPtEtaPhiM(ev.pf_pt[p],ev.pf_eta[p],ev.pf_phi[p],ev.pf_m[p]);
        jet.addParticle(Particle(pp4, ev.pf_c[p], ev.pf_id[p], 0, p, ev.pf_puppiWgt[p]));
        if (ev.pf_c[p] != 0) jet.addTrack(pp4, ev.pf_id[p]);
      }
    }
    
    jets.push_back(jet);
  }
  
  //additional jet-jet information
  for (unsigned int i = 0; i < jets.size(); i++) {
    for (unsigned int j = i+1; j < jets.size(); j++) {
      //flag jet-jet overlaps
      if (jets[i].p4().DeltaR(jets[j].p4()) < 0.8) {
        jets[i].setOverlap(1);
        jets[j].setOverlap(1);
      }
      //flag non-b jets as part of W boson candidates: flavor 0->1
      if (jets[i].flavor()==5 or jets[j].flavor()==5) continue;
      TLorentzVector wCand = jets[i].p4() + jets[j].p4();
      if (abs(wCand.M()-80.4) < 15.) {
        jets[i].setFlavor(1);
        jets[j].setFlavor(1);
      }
    }
  }
  
  return jets;
}

//
std::vector<Jet> SelectionTool::getGenJets(MiniEvent_t ev, double minPt, double maxEta, std::vector<Particle> leptons) {
  std::vector<Jet> jets;
  
  for (int i = 0; i < ev.ng; i++) {
    if (abs(ev.g_id[i])<10) {
      TLorentzVector jp4;
      jp4.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);

      //cross clean with leptons
      bool overlapsWithLepton(false);
      for (auto& lepton : leptons) {
        if(jp4.DeltaR(lepton.p4())<0.4) overlapsWithLepton=true;
      }
      if(overlapsWithLepton) continue;
      
      //jet kinematic selection
      if(jp4.Pt() < minPt || abs(jp4.Eta()) > maxEta) continue;

      //flavor
      int flavor = ev.g_id[i];
      
      Jet jet(jp4, flavor, i);
      
      //fill jet constituents
      for (int p = 0; p < ev.ngpf; p++) {
        if (ev.gpf_g[p] == i) {
          TLorentzVector pp4;
          pp4.SetPtEtaPhiM(ev.gpf_pt[p],ev.gpf_eta[p],ev.gpf_phi[p],ev.gpf_m[p]);
          jet.addParticle(Particle(pp4, ev.gpf_c[p], ev.gpf_id[p], 0, p, 1.));
          if (ev.gpf_c[p] != 0) jet.addTrack(pp4, ev.gpf_id[p]);
        }
      }
      
      jets.push_back(jet);
    }
  }
  
  //additional jet-jet information
  for (unsigned int i = 0; i < jets.size(); i++) {
    for (unsigned int j = i+1; j < jets.size(); j++) {
      //flag jet-jet overlaps
      if (jets[i].p4().DeltaR(jets[j].p4()) < 0.8) {
        jets[i].setOverlap(1);
        jets[j].setOverlap(1);
      }
      //flag non-b jets as part of W boson candidates: flavor 0->1
      if (jets[i].flavor()==5 or jets[j].flavor()==5) continue;
      TLorentzVector wCand = jets[i].p4() + jets[j].p4();
      if (abs(wCand.M()-80.4) < 15.) {
        jets[i].setFlavor(1);
        jets[j].setFlavor(1);
      }
    }
  }
  
  return jets;
}

//
TString SelectionTool::flagFinalState(MiniEvent_t ev, std::vector<Particle> preselleptons) {

  //clear vectors
  leptons_.clear(); 
  vetoLeptons_.clear();
  jets_.clear();

  //if no set of pre-selected leptons has been passed, use standard top selections
  if(preselleptons.size()==0) preselleptons=getTopFlaggedLeptons(ev);

  //check if triggers have fired
  bool hasETrigger(((ev.triggerBits>>0)&0x1)!=0);
  bool hasMTrigger(((ev.triggerBits>>1)&0x3)!=0);
  bool hasEMTrigger(hasETrigger || hasMTrigger || ((ev.triggerBits>>3)&0x8)!=0);
  bool hasMMTrigger(hasETrigger || ((ev.triggerBits>>9)&0x1)!=0);
  bool hasEETrigger(hasETrigger || ((ev.triggerBits>>10)&0x1)!=0);

  //decide the channel based on the lepton multiplicity and set lepton collections
  std::vector<Particle> passllid( getGoodLeptons(preselleptons,PASSLLID) ), passlid( getGoodLeptons(preselleptons,PASSLID) );
  TString chTag("");
  if(passllid.size()>=2){
    int ch( abs(passllid[0].id()*passllid[1].id()) );
    if      (ch==11*13 && hasEMTrigger) chTag = "EM";
    else if (ch==13*13 && hasMMTrigger) chTag = "MM";
    else if (ch==11*11 && hasEETrigger) chTag = "EE";
    leptons_=preselleptons;
  }
  else if(passlid.size()==1){
    int ch(abs(passlid[0].id()) );
    if      (ch==13 && hasMTrigger) chTag = "M";
    else if (ch==11 && hasETrigger) chTag = "E";
    leptons_=preselleptons;
    vetoLeptons_=getGoodLeptons(preselleptons,SelectionTool::PASSLVETO, 0., 99., &leptons_);
  }

  //select jets based on the leptons
  jets_=getGoodJets(ev,30.,2.4,leptons_);
  
  //final consistency check for data 
  if(ev.isData)
    {
      if(chTag=="EE"      && !acceptEE_) chTag="";
      else if(chTag=="MM" && !acceptMM_) chTag="";
      else if(chTag=="EM" && !acceptEM_) chTag="";
      else if(chTag=="E"  && !acceptE_)  chTag="";
      else if(chTag=="M"  && !acceptM_)  chTag="";
    }
 
  //all done
  return chTag;
}

TString SelectionTool::flagGenFinalState(std::vector<Particle> leptons) {
  //decide the channel
  TString chTag("");
  if(leptons.size()>=2) {
    if      (abs(leptons[0].id()*leptons[1].id())==11*13) chTag = "EM";
    else if (abs(leptons[0].id()*leptons[1].id())==13*13) chTag = "MM";
    else if (abs(leptons[0].id()*leptons[1].id())==11*11) chTag = "EE";
  }
  if(leptons.size()==1) {
    if      (abs(leptons[0].id())==13) chTag = "M";
    else if (abs(leptons[0].id())==11) chTag = "E";
  }
  
  return chTag;
}

