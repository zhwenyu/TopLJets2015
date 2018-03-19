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
#include "TopLJets2015/TopAnalysis/interface/VBFVectorBoson.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"


using namespace std;


//
void RunVBFVectorBoson(TString filename,
                       TString outname,
                       Int_t channelSelection, 
                       Int_t chargeSelection, 
                       TH1F *normH, 
                       TString era,
                       Bool_t debug)
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  
  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *genPU=(TH1 *)f->Get("analysis/putrue");
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;

  //LUMINOSITY+PILEUP
  LumiTools lumi(era,genPU);
  
  //LEPTON EFFICIENCIES
  EfficiencyScaleFactorsWrapper gammaEffWR(filename.Contains("Data13TeV"),era);

  //B-TAG CALIBRATION
  BTagSFUtil btvSF(era,"DeepCSV",BTagEntry::OperatingPoint::OP_MEDIUM,"",0);
  
  //JEC/JER
  JECTools jec(era);
  
  //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);
  TString cats[]={"VBFA","HighPtA","MM"};
  ht.addHist("puwgtctr",     new TH1F("puwgtctr",      ";Weight sums;Events",2,0,2));  
  ht.addHist("category",     new TH1F("category",      ";Category;Events",6,0,6));  
  for(size_t i=0; i<3; i++)
    {
      ht.addHist(cats[i]+"_nvtx",         new TH1F(cats[i]+"_nvtx",             ";Vertex multiplicity;Events",100,-0.5,101.5));
      ht.addHist(cats[i]+"_vpt", 	  new TH1F(cats[i]+"_vectorbosonPt",    ";Boson p_{T}[GeV];Events",25,50,500));
      ht.addHist(cats[i]+"_vy", 	  new TH1F(cats[i]+"_vectorbosony",     ";Boson rapidity;Events",25,0,3));
      ht.addHist(cats[i]+"_vystar", 	  new TH1F(cats[i]+"_vectorbosonystar", ";y-(1/2)(y_{j1}+y_{j2});Events",25,0,3));
      ht.addHist(cats[i]+"_njets",        new TH1F(cats[i]+"_njets",            ";Jet multiplicity;Events",15,-0.5,14.5));
      ht.addHist(cats[i]+"_mjj", 	  new TH1F(cats[i]+"_mjj",              ";Dijet invariant mass [GeV];Events",80,400,2000));
      ht.addHist(cats[i]+"_leadpt",       new TH1F(cats[i]+"_leadpt",           ";Leading jet p_{T} [GeV];Events",25,50,500));
      ht.addHist(cats[i]+"_subleadpt",    new TH1F(cats[i]+"_subleadpt"   ,     ";Sub-leading jet p_{T} [GeV];Events",25,50,500));
      ht.addHist(cats[i]+"_centraleta",   new TH1F(cats[i]+"_centraleta",       ";Most central jet |#eta|;Events",25,0,5));
      ht.addHist(cats[i]+"_forwardeta",   new TH1F(cats[i]+"_forwardeta",       ";Most forward jet |#eta|;Events",25,0,5));
      ht.addHist(cats[i]+"_dijetpt",      new TH1F(cats[i]+"_dijetpt",          ";Dijet p_{T} [GeV];Events",25,50,500));
      ht.addHist(cats[i]+"_detajj", 	  new TH1F(cats[i]+"_detajj" ,          ";#Delta#eta(J,J);Events",20,-5,5));
      ht.addHist(cats[i]+"_ht", 	  new TH1F(cats[i]+"_ht",               ";H_{T} [GeV];Events",50,0,500));
    }

  std::cout << "init done" << std::endl;

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, debug, triggerList,SelectionTool::VBF);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);

      
      //assign randomly a run period
      TString period = lumi.assignRunPeriod();
      
      //////////////////
      // CORRECTIONS //
      ////////////////      
      jec.smearJetEnergies(ev);
             
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      TString chTag = selector.flagFinalState(ev);
      if(chTag!="A" && chTag!="MM") continue;
      std::vector<Particle> &photons     = selector.getSelPhotons(); 
      std::vector<Particle> &leptons     = selector.getSelLeptons(); 
      std::vector<Jet>      &jets        = selector.getJets();  

      //boson kinematics
      //for the photon refine also the category according to the bit
      bool isHighPtA(false),isVBFA(false);
      TLorentzVector boson(0,0,0,0);
      if(chTag=="A") {
        isVBFA=(selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v", ev.triggerBits) && photons[0].Pt()>75);
        isHighPtA=( selector.hasTriggerBit("HLT_Photon200_v", ev.triggerBits) && photons[0].Pt()>200 );
        boson += photons[0];
      }
      else {
        boson += leptons[0];
        boson += leptons[1];
      }
      
      //further selection
      bool passBoson(boson.Pt()>75 && fabs(boson.Rapidity())<1.442);
      bool passJets(jets.size()>=2);

      //jet related variables
      float ystar(0),scalarht(0);
      for(auto j : jets) scalarht += j.Pt();
      if(passJets) ystar=boson.Rapidity()-0.5*(jets[0].Rapidity()+jets[1].Rapidity());
      
      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0);
      std::vector<double>plotwgts(1,wgt);
      ht.fill("puwgtctr",0,plotwgts);
      if (!ev.isData) {

        // norm weight
        wgt  = (normH? normH->GetBinContent(1) : 1.0);
        
        // pu weight
        double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
        std::vector<double>puPlotWgts(1,puWgt);
        ht.fill("puwgtctr",1,puPlotWgts);
        
        // photon trigger*selection weights
        float trigSF(1.0), selSF(1.0);
        if(chTag=="A")
          {
            trigSF *= gammaEffWR.getTriggerCorrection({},photons,{}, period).first;
            selSF  *= gammaEffWR.getOfflineCorrection(photons[0], period).first;
          }
        else
          {
            trigSF *=gammaEffWR.getTriggerCorrection(leptons,{},{}, period).first;
            selSF  *=gammaEffWR.getOfflineCorrection(leptons[0], period).first;
            selSF  *=gammaEffWR.getOfflineCorrection(leptons[1], period).first;
          }
        wgt *= puWgt*trigSF*selSF;
        
       
        // generator level weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);

        //update weight for plotter
        plotwgts[0]=wgt;
      }

      //control histograms
      TString c(chTag);
      if(chTag=="A") {


        if(isHighPtA) c="HighPtA";
        if(isVBFA)    c="VBFA";

        ht.fill("category",0,plotwgts);
        if(isHighPtA) ht.fill("category",1,plotwgts);
        if(isHighPtA && isVBFA) ht.fill("category",2,plotwgts);
        if(passBoson && passJets) {
          ht.fill("category",3,plotwgts);
          if(isHighPtA) ht.fill("category",4,plotwgts);
          if(isHighPtA &&isVBFA)ht.fill("category",5,plotwgts);
        }
      }

      ht.fill(c+"_nvtx",  ev.nvtx,          plotwgts);
      ht.fill(c+"_njets", jets.size(),      plotwgts);
      ht.fill(c+"_ht",    scalarht,         plotwgts);
      ht.fill(c+"_vpt",   boson.Pt(),       plotwgts);
      ht.fill(c+"_vy",    boson.Rapidity(), plotwgts);   
      if(passBoson && passJets) {
        ht.fill(c+"_vystar",     ystar, plotwgts);
        ht.fill(c+"_mjj", 	 (jets[0]+jets[1]).M(),	plotwgts);
        ht.fill(c+"_leadpt",     jets[0].Pt(),	plotwgts);
        ht.fill(c+"_subleadpt",  jets[1].Pt(),	plotwgts);
        ht.fill(c+"_centraleta", min(fabs(jets[0].Eta()),fabs(jets[1].Eta())),	plotwgts);
        ht.fill(c+"_forwardeta", max(fabs(jets[0].Eta()),fabs(jets[1].Eta())),	plotwgts);
        ht.fill(c+"_dijetpt",    (jets[0]+jets[1]).Pt(),plotwgts);
        ht.fill(c+"_detajj",     jets[0].Eta()-jets[1].Eta(),plotwgts);
      }
    }
  
  //close input file
  f->Close();
  
  //save histos to file  
  fOut->cd();
  for (auto& it : ht.getPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}
