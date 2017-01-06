#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"


//apply jet energy resolutions
MiniEvent_t smearJetEnergies(MiniEvent_t ev, std::string option) {
  for (int k = 0; k < ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    //smear jet energy resolution for MC
    float genJet_pt(0);
    if(ev.j_g[k]>-1) genJet_pt = ev.g_pt[ ev.j_g[k] ];
    if(!ev.isData && genJet_pt>0) {
      int smearIdx(0);
      if(option=="up") smearIdx=1;
      if(option=="down") smearIdx=2;
      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt)[smearIdx];
      jp4 *= jerSmear;
      ev.j_pt[k]   = jp4.Pt();
      ev.j_eta[k]  = jp4.Eta();
      ev.j_phi[k]  = jp4.Phi();
      ev.j_mass[k] = jp4.M();
    }
  }
  
  return ev;
}

//see working points in https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco
MiniEvent_t addBTagDecisions(MiniEvent_t ev,float wp) {
  for (int k = 0; k < ev.nj; k++) {
    ev.j_btag[k] = (ev.j_csv[k] > wp);
  }
  
  return ev;
}


//details in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
MiniEvent_t updateBTagDecisions(MiniEvent_t ev, 
				std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> &btvsfReaders,
				std::map<BTagEntry::JetFlavor, TGraphAsymmErrors*> &expBtagEff, 
				std::map<BTagEntry::JetFlavor, TGraphAsymmErrors*> &expBtagEffPy8, 
				BTagSFUtil *myBTagSFUtil, 
				std::string option) {
  for (int k = 0; k < ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    bool isBTagged(ev.j_btag[k]);
    if(!ev.isData) {
      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
      float expEff(1.0), jetBtagSF(1.0);
      
      BTagEntry::JetFlavor hadFlav=BTagEntry::FLAV_UDSG;
      if(abs(ev.j_hadflav[k])==4) hadFlav=BTagEntry::FLAV_C;
      if(abs(ev.j_hadflav[k])==5) hadFlav=BTagEntry::FLAV_B;

      expEff    = expBtagEff[hadFlav]->Eval(jptForBtag); 
      jetBtagSF = btvsfReaders[hadFlav]->eval_auto_bounds( option, hadFlav, jetaForBtag, jptForBtag);
      jetBtagSF *= expEff>0 ? expBtagEffPy8[hadFlav]->Eval(jptForBtag)/expBtagEff[hadFlav]->Eval(jptForBtag) : 0.;
      
      //updated b-tagging decision with the data/MC scale factor
      myBTagSFUtil->modifyBTagsWithSF(isBTagged, jetBtagSF, expEff);
      ev.j_btag[k] = isBTagged;
    }
  }
  
  return ev;
}

//details in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> getBTVcalibrationReaders(TString era,BTagEntry::OperatingPoint btagOP)
{
  //start the btag calibration
  TString btagUncUrl(era+"/btagSFactors.csv");
  gSystem->ExpandPathName(btagUncUrl);
  BTagCalibration btvcalib("csvv2", btagUncUrl.Data());

  //start calibration readers for b,c and udsg separately including the up/down variations
  std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> btvCalibReaders;
  btvCalibReaders[BTagEntry::FLAV_B]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders[BTagEntry::FLAV_B]->load(btvcalib,BTagEntry::FLAV_B,"mujets");
  btvCalibReaders[BTagEntry::FLAV_C]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders[BTagEntry::FLAV_C]->load(btvcalib,BTagEntry::FLAV_C,"mujets");
  btvCalibReaders[BTagEntry::FLAV_UDSG]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders[BTagEntry::FLAV_UDSG]->load(btvcalib,BTagEntry::FLAV_UDSG,"incl");

  //all done
  return btvCalibReaders;
}

//the expections are created with the script scripts/saveExpectedBtagEff.py (cf README)
std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> readExpectedBtagEff(TString era,TString btagExpPostFix)
{
  //open up the ROOT file with the expected efficiencies
  TString btagEffExpUrl(era+"/expTageff.root");
  btagEffExpUrl.ReplaceAll(".root",btagExpPostFix+".root");
  gSystem->ExpandPathName(btagEffExpUrl);
  TFile *beffIn=TFile::Open(btagEffExpUrl);
  
  //read efficiency graphs
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> expBtagEff;
  expBtagEff[BTagEntry::FLAV_B]=(TGraphAsymmErrors *)beffIn->Get("b");
  expBtagEff[BTagEntry::FLAV_C]=(TGraphAsymmErrors *)beffIn->Get("c");
  expBtagEff[BTagEntry::FLAV_UDSG]=(TGraphAsymmErrors *)beffIn->Get("udsg");
  beffIn->Close();

  //all done
  return expBtagEff;
}

