#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-17-010.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include "TMath.h"

using namespace std;

//
void TOP17010::init(){

  //read input
  f_           = TFile::Open(filename_);
  isSignal_    = filename_.Contains("_TT");
  triggerList_ = (TH1F *)f_->Get("analysis/triggerList");
  weightSysts_ = getWeightSysts(f_);
  t_           = (TTree*)f_->Get("analysis/data");
  attachToMiniEventTree(t_,ev_,true);
  nentries_   = t_->GetEntriesFast();
  if (debug_) nentries_ = 10000; //restrict number of entries for testing
  t_->GetEntry(0);


  baseName_=gSystem->BaseName(outname_); 
  TString dirName=gSystem->DirName(outname_);

  //corrections
  lumi_ = new LumiTools(era_,genPU_);

  std::map<TString,TString> cfgMap;
  cfgMap["g_id"]     = "Tight";
  cfgMap["m_id"]     = "TightID";
  cfgMap["m_iso"]    = "TightRelIso";
  cfgMap["m_id4iso"] = "TightIDandIPCut";
  cfgMap["e_id"]     = "MVA80";
  gammaEffWR_  = new EfficiencyScaleFactorsWrapper(filename_.Contains("Data13TeV"),era_,cfgMap);
  l1PrefireWR_ = new L1PrefireEfficiencyWrapper(filename_.Contains("Data13TeV"),era_);  
  btvSF_       = new BTagSFUtil(era_);

  //theory uncs
  fragWeights_   = getBFragmentationWeights(era_);
  semilepBRwgts_ = getSemilepBRWeights(era_);

  //start selection tool
  selector_ = new SelectionTool(filename_, debug_, triggerList_);
}


//
void TOP17010::bookHistograms() {

  ht_ = new HistTool(0);
  ht_->addHist("puwgtctr",   new TH1F("puwgtctr", ";Weight sums;Events",                        2,0,2));  
  ht_->addHist("nvtx",       new TH1F("nvtx",     ";Vertex multiplicity;Events",                100,-0.5,99.5));
  ht_->addHist("mll", 	     new TH1F("mll",      ";Dilepton invariant mass [GeV];Events",      50,0,550));  
  ht_->addHist("ptll", 	     new TH1F("ptll",     ";Dilepton p_{T}[GeV];Events",                50,0,550));  
  ht_->addHist("l1pt",       new TH1F("l1pt",     ";Lepton 1 transverse momentum [GeV];Events", 50,20,200));
  ht_->addHist("l1eta",      new TH1F("l1eta",    ";Lepton 1 pseudo-rapidity;Events",           10,0,2.5));
  ht_->addHist("l2pt",       new TH1F("l2pt",     ";Lepton 2 transverse momentum [GeV];Events", 50,20,200));
  ht_->addHist("l2eta",      new TH1F("l2eta",     ";Lepton 2 pseudo-rapidity;Events",          10,0,2.5));

  ht_->addHist("njets",      new TH1F("njets",    ";Jet multiplicity;Events",                 6,0,6)); 
  ht_->addHist("nbjets",     new TH1F("nbjets",   ";b jet multiplicity;Events",               5,0,5));
  ht_->addHist("j1pt",       new TH1F("j1pt",     ";Jet 1 transverse momentum [GeV];Events",  50,30,200));
  ht_->addHist("j1eta",      new TH1F("j1eta",    ";Jet 1 pseudo-rapidity;Events",            10,0,2.5));
  ht_->addHist("j2pt",       new TH1F("j2pt",     ";Jet 2 transverse momentum [GeV];Events",  50,30,200));
  ht_->addHist("j2eta",      new TH1F("j2eta",    ";Jet 2 pseudo-rapidity;Events",            10,0,2.5));

  ht_->addHist("evcount",    new TH1F("evcount",  ";Pass;Events", 1,0,1));  
  ht_->addHist("drlb",       new TH1F("drlb",     ";#DeltaR(lb);Events",      50,0,2*TMath::Pi()) );

  TFile *rIn=TFile::Open("$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/top17010/mlbresol.root");  
  std::vector<TString> templates={"mlb","ptlb"};
  for(auto t : templates) {
    TH1 *th=(TH1 *)rIn->Get(t);
    th->SetDirectory(0);
    th->Reset("ICE");
    ht_->addHist(t,th);  
  }
  rIn->Close();

  TString expSystNames[]={"puup",        "pudn",
                          "eetrigup",    "eetrigdn",
                          "eselup",      "eseldn",
                          "mselup",      "mseldn",
                          "l1prefireup", "l1prefiredn",
                          "ees1up", "ees1dn", "ees2up", "ees2dn", "ees3up", "ees3dn", "ees4up", "ees4dn",  "ees5up", "ees5dn",  "ees6up", "ees6dn",  "ees7up", "ees7dn",
                          "mes1up", "mes1dn", "mes2up", "mes2dn", "mes3up", "mes3dn", "mes4up", "mes4dn"
                          "JERup",       "JERdn",
                          "topptup",     "topptdn",
                          "AbsoluteStatJECup","AbsoluteScaleJECup","AbsoluteMPFBiasJECup","FragmentationJECup","SinglePionECALJECup","SinglePionHCALJECup","FlavorPureGluonJECup","FlavorPureQuarkJECup","FlavorPureCharmJECup","FlavorPureBottomJECup","TimePtEtaJECup","RelativeJEREC1JECup","RelativeJEREC2JECup","RelativeJERHFJECup","RelativePtBBJECup","RelativePtEC1JECup","RelativePtEC2JECup","RelativePtHFJECup","RelativeBalJECup","RelativeFSRJECup","RelativeStatFSRJECup","RelativeStatECJECup","RelativeStatHFJECup","PileUpDataMCJECup","PileUpPtRefJECup","PileUpPtBBJECup","PileUpPtEC1JECup","PileUpPtEC2JECup","PileUpPtHFJECup",
                          "AbsoluteStatJECdn","AbsoluteScaleJECdn","AbsoluteMPFBiasJECdn","FragmentationJECdn","SinglePionECALJECdn","SinglePionHCALJECdn","FlavorPureGluonJECdn","FlavorPureQuarkJECdn","FlavorPureCharmJECdn","FlavorPureBottomJECdn","TimePtEtaJECdn","RelativeJEREC1JECdn","RelativeJEREC2JECdn","RelativeJERHFJECdn","RelativePtBBJECdn","RelativePtEC1JECdn","RelativePtEC2JECdn","RelativePtHFJECdn","RelativeBalJECdn","RelativeFSRJECdn","RelativeStatFSRJECdn","RelativeStatECJECdn","RelativeStatHFJECdn","PileUpDataMCJECdn","PileUpPtRefJECdn","PileUpPtBBJECdn","PileUpPtEC1JECdn","PileUpPtEC2JECdn","PileUpPtHFJECdn",
                          "bfragup", "bfragdn"
                          "slbrup",  "slbrdn"};
  
  //instantiate 2D histograms for most relevant variables to trace with systs
  TString hoi[]={"ptlb","drlb","mlb","evcount"};
  size_t nexpSysts=sizeof(expSystNames)/sizeof(TString);
  expSysts_=std::vector<TString>(expSystNames,expSystNames+nexpSysts);  
  for(size_t ih=0; ih<sizeof(hoi)/sizeof(TString); ih++)
    {
      TH1 *histo=ht_->getPlots()[hoi[ih]];
      
      //experimental systs
      ht_->addHist(hoi[ih]+"_exp",      
                  new TH2F(hoi[ih]+"_exp", 
                           Form(";%s;Experimental systematic variation;Events",histo->GetName()),
                           histo->GetNbinsX(),histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax(),
                           nexpSysts,0,nexpSysts));
      for(size_t is=0; is<nexpSysts; is++)
        ht_->get2dPlots()[hoi[ih]+"_exp"]->GetYaxis()->SetBinLabel(is+1,expSystNames[is]);
      
      //theory systs
      size_t nthSysts(weightSysts_.size()); 
      if(nthSysts>0){
        ht_->addHist(hoi[ih]+"_th",      
                    new TH2F(hoi[ih]+"_th", 
                             Form(";%s;Theory systematic variation;Events",histo->GetName()),
                             histo->GetNbinsX(),histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax(),
                             nthSysts,0,nthSysts));
        for(size_t is=0; is<nthSysts; is++)
          ht_->get2dPlots()[hoi[ih]+"_th"]->GetYaxis()->SetBinLabel(is+1,weightSysts_[is].first);
      }
    }
}



//
void TOP17010::runAnalysis()
{
  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  for (Int_t iev=0;iev<nentries_;iev++)
    {
      t_->GetEntry(iev);
      if(debug_) cout << "Number of event: "<<iev<<endl;
      if(iev%10000==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries_);
       
      //////////////////
      // CORRECTIONS  //
      //////////////////
      TString period = lumi_->assignRunPeriod();
      btvSF_->addBTagDecisions(ev_);
      if(!ev_.isData) btvSF_->updateBTagDecisions(ev_);      
      
      //TRIGGER
      bool hasETrigger(  selector_->hasTriggerBit("HLT_Ele32_eta2p1_WPTight_Gsf_v",                       ev_.triggerBits) );
      bool hasMTrigger(  selector_->hasTriggerBit("HLT_IsoMu24_v",                                        ev_.triggerBits) || 
                         selector_->hasTriggerBit("HLT_IsoTkMu24_v",                                      ev_.triggerBits) );
      bool hasEMTrigger( selector_->hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     ev_.triggerBits) ||
                         selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",    ev_.triggerBits) ||
                         selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits) ||
                         selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev_.triggerBits) ||
                         selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits) );
      bool hasMMTrigger( selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",                ev_.triggerBits) ||
                         selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",              ev_.triggerBits) );
      bool hasEETrigger( selector_->hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          ev_.triggerBits) );
      if(era_.Contains("2017")) {
        hasETrigger =(selector_->hasTriggerBit("HLT_Ele35_WPTight_Gsf_v", ev_.triggerBits));
        hasMTrigger =(selector_->hasTriggerBit("HLT_IsoMu24_v",           ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_IsoMu24_2p1_v",       ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_IsoMu27_v",           ev_.triggerBits) );          
        hasEMTrigger=(selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",    ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev_.triggerBits) );
        hasMMTrigger=(selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                  ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",          ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",        ev_.triggerBits) );
        hasEETrigger=(selector_->hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",             ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          ev_.triggerBits) );
      }
      if(ev_.isData) {
        if(era_.Contains("2016")) {
          hasEMTrigger=false;
          if(ev_.run<=280385)
            {
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev_.triggerBits);
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",  ev_.triggerBits);
            }
          if(ev_.run>=278273 && ev_.run<=280385)
            {
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", ev_.triggerBits);
            }
          if(ev_.run>=278273)
            {
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev_.triggerBits);
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits);
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits);
            }
        }
        if(filename_.Contains("2017E") || filename_.Contains("2017F")){
          hasMTrigger=selector_->hasTriggerBit("HLT_IsoMu27_v",ev_.triggerBits);
        }
        if(!(filename_.Contains("2017A") || filename_.Contains("2017B"))){
          hasMMTrigger=(selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",   ev_.triggerBits) ||
                        selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev_.triggerBits) );
        }
        
        //trigger consistency with PD
        if( selector_->isDoubleEGPD()       && !hasEETrigger ) return;
        if( selector_->isDoubleMuonPD()     && !hasMMTrigger ) return;
        if( selector_->isMuonEGPD()         && !hasEMTrigger ) return;
        if( selector_->isSingleElectronPD() && !hasETrigger )  return;
        if( selector_->isSingleMuonPD()     && !hasMTrigger )  return;          
      }
    
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      std::vector<Particle> leptons = selector_->flaggedLeptons(ev_);     
      leptons = selector_->selLeptons(leptons,SelectionTool::MEDIUM,SelectionTool::MVA80,20,2.5);
      std::vector<Jet> alljets = selector_->getGoodJets(ev_,30.,4.7,leptons);
      TopWidthEvent twe(leptons,alljets);
      if(twe.dilcode==0) {
        return;
      }
      else if(twe.dilcode==11*11) {
        if(!(hasEETrigger || hasETrigger)) return;
        if(selector_->isSingleElectronPD() && hasEETrigger) return;
      }
      else if(twe.dilcode==13*13) {
        if(!(hasMMTrigger || hasMTrigger)) return;
        if(selector_->isSingleMuonPD() && hasMMTrigger) return;
      }
      else if(twe.dilcode==11*13) {
        if(!(hasEMTrigger || hasETrigger || hasMTrigger)) return;
        if(selector_->isSingleMuonPD()     && hasEMTrigger) return;
        if(selector_->isSingleElectronPD() && (hasEMTrigger || hasMTrigger)) return;
      }
       
      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0);
      std::vector<float>puWgts(3,1.0),topptWgts(2,1.0),bfragWgts(2,1.0),slbrWgts(2,1.0);
      EffCorrection_t trigSF(1.0,0.),l1SF(1.0,0.),l2SF(1.0,0.0),l1trigprefireProb(1.0,0.);
      if (!ev_.isData) {
            
        // pu weight
        std::vector<double> plotwgts(1,1);
        ht_->fill("puwgtctr",0,plotwgts);
        puWgts=lumi_->pileupWeight(ev_.g_pu,period);
        std::vector<double>puPlotWgts(1,puWgts[0]);
        ht_->fill("puwgtctr",1,puPlotWgts);

        //L1 prefire probability
        //        l1trigprefireProb=l1PrefireWR_->getCorrection(alljets);

        // photon trigger*selection weights        
        trigSF = gammaEffWR_->getTriggerCorrection(leptons,{},{}, period);
        l1SF   = gammaEffWR_->getOfflineCorrection(leptons[0].id(),leptons[0].pt(),leptons[0].eta(), period);
        l2SF   = gammaEffWR_->getOfflineCorrection(leptons[1].id(),leptons[1].pt(),leptons[1].eta(), period);

        //for signal top pt weights        
        if(isSignal_) {
          for(Int_t igen=0; igen<ev_.ngtop; igen++)
            {
              if(abs(ev_.gtop_id[igen])!=6) continue;
              float topsf=TMath::Exp(0.156-0.00137*ev_.gtop_pt[igen]);
              topptWgts[0] *= topsf;
              topptWgts[1] *= 1./topsf;
            }
        }

        //b-fragmentation and semi-leptonic branching fractions
        bfragWgts[0] = computeBFragmentationWeight(ev_,fragWeights_["downFrag"]);
        bfragWgts[1] = computeBFragmentationWeight(ev_,fragWeights_["upFrag"]);
        slbrWgts[0]  = computeSemilepBRWeight(ev_,semilepBRwgts_["semilepbrDown"]);
        slbrWgts[1]  = computeSemilepBRWeight(ev_,semilepBRwgts_["semilepbrUp"]);        

        //final nominal weight
        wgt *= (normH_? normH_->GetBinContent(1) : 1.0);
        wgt *= (ev_.g_nw>0 ? ev_.g_w[0] : 1.0);
        wgt *= puWgts[0]*l1trigprefireProb.first*trigSF.first*l1SF.first*l2SF.second;
      }  
      fillControlHistograms(twe,wgt);

      //FIXME

      //experimental systs cycle: better not to do anything else after this...
      //final category selection is repeated ad nauseam with varied objects/weights and mva is re-evaluated several times
      if(ev_.isData) continue;
      selector_->setDebug(false);
      for(size_t is=0; is<expSysts_.size(); is++){

        //uncertainty
        TString sname=expSysts_[is];
        bool isUpVar(sname.Contains("up"));
        
        //base values and kinematics
        TString icat(twe.cat);
        float iwgt=(ev_.g_nw>0 ? ev_.g_w[0] : 1.0);
        iwgt *= (normH_? normH_->GetBinContent(1) : 1.0);

        bool reSelect(false);        

        //FIXME
        EffCorrection_t selSF(1.0,0.0);
        if(sname=="puup")             iwgt *= puWgts[1]*trigSF.first*selSF.first*l1trigprefireProb.first;
        else if(sname=="pudn")        iwgt *= puWgts[2]*trigSF.first*selSF.first*l1trigprefireProb.first;
        else if(sname=="trigup")      iwgt *= puWgts[0]*(trigSF.first+trigSF.second)*selSF.first*l1trigprefireProb.first;
        else if(sname=="trigdn")      iwgt *= puWgts[0]*(trigSF.first-trigSF.second)*selSF.first*l1trigprefireProb.first;
        else if(sname=="selup")       iwgt *= puWgts[0]*trigSF.first*(selSF.first+selSF.second)*l1trigprefireProb.first;
        else if(sname=="seldn")       iwgt *= puWgts[0]*trigSF.first*(selSF.first-selSF.second)*l1trigprefireProb.first;
        else if(sname=="l1prefireup") iwgt *= puWgts[0]*trigSF.first*selSF.first*(l1trigprefireProb.first+l1trigprefireProb.second);
        else if(sname=="l1prefiredn") iwgt *= puWgts[0]*trigSF.first*selSF.first*(l1trigprefireProb.first-l1trigprefireProb.second);
        else if(sname=="topptup")     iwgt = wgt*topptWgts[0];
        else if(sname=="topptdn")     iwgt = wgt*topptWgts[1];
        else if(sname=="bfragup")     iwgt = wgt*bfragWgts[0];
        else if(sname=="bfragdn")     iwgt = wgt*bfragWgts[1];
        else if(sname=="slbrup")      iwgt = wgt*slbrWgts[0];
        else if(sname=="slbrdn")      iwgt = wgt*slbrWgts[1];
        else                          iwgt = wgt;           

        std::vector<Jet> ijets(alljets);

        //FIXME
        if(sname.Contains("JEC") || sname.Contains("JER") )  {
          reSelect=true;
          int jecIdx=-1;
          if(sname.Contains("AbsoluteStat"))     jecIdx=0;
          if(sname.Contains("AbsoluteScale"))    jecIdx=1; 
          if(sname.Contains("AbsoluteMPFBias"))  jecIdx=2; 
          if(sname.Contains("Fragmentation"))    jecIdx=3; 
          if(sname.Contains("SinglePionECAL"))   jecIdx=4; 
          if(sname.Contains("SinglePionHCAL"))   jecIdx=5; 
          if(sname.Contains("FlavorPureGluon"))  jecIdx=6; 
          if(sname.Contains("FlavorPureQuark"))  jecIdx=7; 
          if(sname.Contains("FlavorPureCharm"))  jecIdx=8; 
          if(sname.Contains("FlavorPureBottom")) jecIdx=9; 
          if(sname.Contains("TimePtEta"))        jecIdx=10; 
          if(sname.Contains("RelativeJEREC1"))   jecIdx=11; 
          if(sname.Contains("RelativeJEREC2"))   jecIdx=12; 
          if(sname.Contains("RelativeJERHF"))    jecIdx=13; 
          if(sname.Contains("RelativePtBB"))     jecIdx=14; 
          if(sname.Contains("RelativePtEC1"))    jecIdx=15; 
          if(sname.Contains("RelativePtEC2"))    jecIdx=16; 
          if(sname.Contains("RelativePtHF"))     jecIdx=17; 
          if(sname.Contains("RelativeBal"))      jecIdx=18; 
          if(sname.Contains("RelativeFSR"))      jecIdx=19; 
          if(sname.Contains("RelativeStatFSR"))  jecIdx=20; 
          if(sname.Contains("RelativeStatEC"))   jecIdx=21; 
          if(sname.Contains("RelativeStatHF"))   jecIdx=22; 
          if(sname.Contains("PileUpDataMC"))     jecIdx=23; 
          if(sname.Contains("PileUpPtRef"))      jecIdx=24; 
          if(sname.Contains("PileUpPtBB"))       jecIdx=25; 
          if(sname.Contains("PileUpPtEC1"))      jecIdx=26; 
          if(sname.Contains("PileUpPtEC2"))      jecIdx=27; 
          if(sname.Contains("PileUpPtHF"))       jecIdx=28;
          
          //re-scale and re-select jets
          std::vector<Jet> newJets = selector_->getGoodJets(ev_,20.,2.4,leptons);
          ijets.clear();
          for(auto j : alljets) {

            int idx=j.getJetIndex();

            //shift jet energy
            float scaleVar(1.0);
            if(jecIdx<0) {
              scaleVar=isUpVar ? ev_.j_jerUp[idx] : ev_.j_jerDn[idx];
            } 
            else {
              int jflav(abs(ev_.j_flav[idx]));
              bool flavorMatches(true);
              if(jecIdx==6 && jflav!=21) flavorMatches=false; //FlavorPureGluon
              if(jecIdx==7 && jflav>=4)  flavorMatches=false; //FlavorPureQuark
              if(jecIdx==8 && jflav!=4)  flavorMatches=false; //FlavorPureCharm
              if(jecIdx==9 && jflav!=5)  flavorMatches=false; //FlavorPureGluon
              if(flavorMatches)
                scaleVar=isUpVar ? ev_.j_jecUp[jecIdx][idx] : ev_.j_jecDn[jecIdx][idx];
            }
            j*=scaleVar;
            //FIXME hardcoded
            if(j.Pt()<30. || fabs(j.Eta())>2.4) continue;

            ijets.push_back(j);
          }
        }

        //re-select if needed
        if(reSelect) {
          
          if (ijets.size()<2) continue;
        }

        //fill with new values/weights
        std::vector<double> eweights(1,iwgt);
        //FIXME
      }
      selector_->setDebug(debug_);
    }
  
  //close input file
  f_->Close();


  //save histograms to the output
  fOut_->cd();
  for (auto& it : ht_->getPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut_); 
    it.second->Write(); 
  }
  for (auto& it : ht_->get2dPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut_); 
    it.second->Write(); 
  }
  fOut_->Close();		
}

//
void TOP17010::fillControlHistograms(TopWidthEvent &twe,float &wgt) {

  //plot weight
  std::vector<double> cplotwgts(1,wgt);
  ht_->fill("nvtx",    ev_.nvtx,        cplotwgts,twe.cat);        
  ht_->fill("mll",     twe.mll,         cplotwgts,twe.cat);
  ht_->fill("ptll",    twe.ptll,        cplotwgts,twe.cat);
  ht_->fill("l1pt",    twe.l1pt,        cplotwgts,twe.cat);
  ht_->fill("l1eta",   fabs(twe.l1eta), cplotwgts,twe.cat);
  ht_->fill("l2pt",    twe.l2pt,        cplotwgts,twe.cat);
  ht_->fill("l2eta",   fabs(twe.l2eta), cplotwgts,twe.cat);
  ht_->fill("njets",   twe.njets,       cplotwgts,twe.cat);
  ht_->fill("nbjets",  twe.nbjets,      cplotwgts,twe.cat);
  ht_->fill("j1pt",    twe.j1pt,        cplotwgts,twe.cat);
  ht_->fill("j1eta",   fabs(twe.j2eta), cplotwgts,twe.cat);
  ht_->fill("j2pt",    twe.j2pt,        cplotwgts,twe.cat);
  ht_->fill("j2eta",   fabs(twe.j2eta), cplotwgts,twe.cat);
  ht_->fill("evcount", 0,               cplotwgts,twe.cat);
  std::vector<TLorentzVector> lbPairs=twe.getPairs();
  for(auto p4 : lbPairs) {
    ht_->fill("ptlb", p4.Pt(), cplotwgts, twe.cat);
    ht_->fill("mlb",  p4.M(),  cplotwgts, twe.cat);
    TString c(twe.cat+(p4.Pt()>100?"highpt":"lowpt"));
    ht_->fill("mlb",  p4.M(),  cplotwgts, c);
  }

  //theory uncertainties are filled only for MC
  if(ev_.isData) return;
  if(ev_.g_w[0]==0 || normH_==NULL || normH_->GetBinContent(1)==0) return;
    
  //replicas for theory systs
  for(size_t is=0; is<weightSysts_.size(); is++){
    std::vector<double> sweights(1,cplotwgts[0]);
    size_t idx=weightSysts_[is].second;
    sweights[0] *= (ev_.g_w[idx]/ev_.g_w[0])*(normH_->GetBinContent(idx+1)/normH_->GetBinContent(1));
    for(auto p4 : lbPairs) {
      ht_->fill2D("ptlb_th",  p4.Pt(), is, sweights, twe.cat);
      ht_->fill2D("mlb_th",   p4.M(),  is, sweights, twe.cat);
      TString c(twe.cat+(p4.Pt()>100?"highpt":"lowpt"));
      ht_->fill2D("mlb_th",   p4.M(),  is, sweights, c);
    }
  }
}
