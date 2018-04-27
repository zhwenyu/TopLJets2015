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

using namespace std;

//
void RunVBFVectorBoson(TString filename,
                       TString outname,
                       Int_t channelSelection, 
                       Int_t chargeSelection, 
                       TH1F *normH,
                       TH1F *genPU,
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
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  bool isQCDEMEnriched( filename.Contains("MC13TeV_QCDEM") );

  cout << "...producing " << outname << " from " << nentries << " events" << endl;

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
  ht.addHist("puwgtctr", new TH1F("puwgtctr", ";Weight sums;Events",2,0,2));  
  ht.addHist("category", new TH1F("category", ";Category;Events",3,0,3));  
  ht.addHist("qscale",   new TH1F("qscale",   ";Q^{2} scale;Events",100,0,2000));  
  TString cats[]={"VBFA","HighPtA","VBFMM","HighPtMM"};
  for(size_t i=0; i<4; i++)
    {
      ht.addHist(cats[i]+"_nvtx",         new TH1F(cats[i]+"_nvtx",             ";Vertex multiplicity;Events",100,-0.5,99.5));
      ht.addHist(cats[i]+"_vpt", 	  new TH1F(cats[i]+"_vectorbosonPt",    ";Boson p_{T}[GeV];Events",25,50,500));
      ht.addHist(cats[i]+"_vy", 	  new TH1F(cats[i]+"_vectorbosony",     ";Boson rapidity;Events",25,0,3));
      ht.addHist(cats[i]+"_vystar", 	  new TH1F(cats[i]+"_vectorbosonystar", ";y-(1/2)(y_{j1}+y_{j2});Events",25,0,3));
      ht.addHist(cats[i]+"_njets",        new TH1F(cats[i]+"_njets",            ";Jet multiplicity;Events",15,-0.5,14.5));
      ht.addHist(cats[i]+"_mjj", 	  new TH1F(cats[i]+"_mjj",              ";Dijet invariant mass [GeV];Events",40,400,2400));
      ht.addHist(cats[i]+"_leadpt",       new TH1F(cats[i]+"_leadpt",           ";Leading jet p_{T} [GeV];Events",25,50,500));
      ht.addHist(cats[i]+"_subleadpt",    new TH1F(cats[i]+"_subleadpt"   ,     ";Sub-leading jet p_{T} [GeV];Events",25,50,500));
      ht.addHist(cats[i]+"_centraleta",   new TH1F(cats[i]+"_centraleta",       ";Most central jet |#eta|;Events",25,0,5));
      ht.addHist(cats[i]+"_forwardeta",   new TH1F(cats[i]+"_forwardeta",       ";Most forward jet |#eta|;Events",25,0,5));
      ht.addHist(cats[i]+"_dijetpt",      new TH1F(cats[i]+"_dijetpt",          ";Dijet p_{T} [GeV];Events",25,50,500));
      ht.addHist(cats[i]+"_detajj", 	  new TH1F(cats[i]+"_detajj" ,          ";#Delta#eta(J,J);Events",20,0,6));
      ht.addHist(cats[i]+"_ht", 	  new TH1F(cats[i]+"_ht",               ";H_{T} [GeV];Events",50,0,8300));
      ht.addHist(cats[i]+"_balance", 	  new TH1F(cats[i]+"_balance",          ";System p_{T} balance [GeV];Events",25,0,250));
	  
      // Study of jet variables
	  ht.addHist(cats[i]+"_jet_c1_00", 	  new TH1F(cats[i]+"_jet_c1_00",          ";Jet shape var. c1_00;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c1_02", 	  new TH1F(cats[i]+"_jet_c1_02",          ";Jet shape var. c1_02;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c1_05", 	  new TH1F(cats[i]+"_jet_c1_05",          ";Jet shape var. c1_05;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c1_10", 	  new TH1F(cats[i]+"_jet_c1_10",          ";Jet shape var. c1_10;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c1_20", 	  new TH1F(cats[i]+"_jet_c1_20",          ";Jet shape var. c1_20;Jets",100,-1,1));

	  ht.addHist(cats[i]+"_jet_c2_00", 	  new TH1F(cats[i]+"_jet_c2_00",          ";Jet shape var. c2_00;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c2_02", 	  new TH1F(cats[i]+"_jet_c2_02",          ";Jet shape var. c2_02;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c2_05", 	  new TH1F(cats[i]+"_jet_c2_05",          ";Jet shape var. c2_05;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c2_10", 	  new TH1F(cats[i]+"_jet_c2_10",          ";Jet shape var. c2_10;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c2_20", 	  new TH1F(cats[i]+"_jet_c2_20",          ";Jet shape var. c2_20;Jets",100,-1,1));

	  ht.addHist(cats[i]+"_jet_c3_00", 	  new TH1F(cats[i]+"_jet_c3_00",          ";Jet shape var. c3_00;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c3_02", 	  new TH1F(cats[i]+"_jet_c3_02",          ";Jet shape var. c3_02;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c3_05", 	  new TH1F(cats[i]+"_jet_c3_05",          ";Jet shape var. c3_05;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c3_10", 	  new TH1F(cats[i]+"_jet_c3_10",          ";Jet shape var. c3_10;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_c3_20", 	  new TH1F(cats[i]+"_jet_c3_20",          ";Jet shape var. c3_20;Jets",100,-1,1));

	  ht.addHist(cats[i]+"_jet_zg", 	  new TH1F(cats[i]+"_jet_zg",          ";Jet shape var. zg;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_gaptd", 	  new TH1F(cats[i]+"_jet_gaptd",          ";Jet shape var. gaptd;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_gawidth", 	  new TH1F(cats[i]+"_jet_gawidth",          ";Jet shape var. gawidth;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_gathrust", 	  new TH1F(cats[i]+"_jet_gathrust",          ";Jet shape var. gathrust;Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_tau21", 	  new TH1F(cats[i]+"_jet_tau21",          ";Jet shape var. #tau_{21};Jets",100,-1,1));
	  ht.addHist(cats[i]+"_jet_tau32", 	  new TH1F(cats[i]+"_jet_tau32",          ";Jet shape var. #tau_{32};Jets",100,-1,1));

      //additional variables from https://link.springer.com/content/pdf/10.1140/epjc/s10052-017-5315-6.pdf
      ht.addHist(cats[i]+"_jjetas", 	  new TH1F(cats[i]+"_jjetas",          ";#eta_{j1}#eta_{j2};Events",200,-25,25));
      ht.addHist(cats[i]+"_centjy", 	  new TH1F(cats[i]+"_centjy",          ";Central jet rapidity;Jets",25,0,3));
      ht.addHist(cats[i]+"_ncentj", 	  new TH1F(cats[i]+"_ncentjj",          ";Number of central jets;Events",10,-0.5,9.5));
      ht.addHist(cats[i]+"_dphivj0", 	  new TH1F(cats[i]+"__dphivj0",          ";#Delta#phi(V,j0);Jets",20,0,4));
      ht.addHist(cats[i]+"_dphivj1", 	  new TH1F(cats[i]+"__dphivj1",          ";#Delta#phi(V,j1);Jets",20,0,4));
      ht.addHist(cats[i]+"_dphivj2", 	  new TH1F(cats[i]+"__dphivj2",          ";#Delta#phi(V,j2);Jets",20,0,4));
      ht.addHist(cats[i]+"_dphivj3", 	  new TH1F(cats[i]+"__dphivj3",          ";#Delta#phi(V,j3);Jets",20,0,4));
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
      if(iev%10000==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
	  cout<<"Event Number "<<iev<<endl;
      std::vector<double>plotwgts(1,1.0);
      ht.fill("qscale",ev.g_qscale,plotwgts);
      
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

      //jet related variables and selection
      float mjj(jets.size()>=2 ?  (jets[0]+jets[1]).M() : 0.);
      float scalarht(0.);
      for(auto j : jets) scalarht += j.Pt();
      bool passJets(jets.size()>=2 && mjj>400);

      //categorize the event according to the boson kinematics
      //for the photon refine also the category according to the trigger  bit
      TLorentzVector boson(0,0,0,0);
      bool isHighPt(false),isVBF(false);
      if(chTag=="A") {
        
        boson += photons[0];
        isVBF    = (selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v", ev.triggerBits) 
                    && photons[0].Pt()>75 
                    && fabs(photons[0].Eta())<1.442
                    && mjj>400);
        isHighPt = ( selector.hasTriggerBit("HLT_Photon200_v", ev.triggerBits) 
                     && photons[0].Pt()>200 );       

        //veto prompt photons on the QCDEM enriched sample
        if( isQCDEMEnriched && ev.gamma_isPromptFinalState[ photons[0].originalReference() ] ) {
          isVBF=false;
          isHighPt=false;
        }
          
      }
      else {
        boson += leptons[0];
        boson += leptons[1];
        isVBF    = boson.Pt()>75 && fabs(boson.Rapidity())<1.442 && mjj>400;
        isHighPt = boson.Pt()>200;
      }
      if(!isVBF && !isHighPt) continue;      
      std::vector<TString> chTags;      
      if(isVBF)    chTags.push_back("VBF"+chTag);
      if(isHighPt) chTags.push_back("HighPt"+chTag);
      
      //system variables
      float ystar(0),balance(0);
      if(passJets) {
        ystar=boson.Rapidity()-0.5*(jets[0].Rapidity()+jets[1].Rapidity());
        balance=(boson+jets[0]+jets[1]).Pt();
      }

      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0);
      if (!ev.isData) {

        // norm weight
        wgt  = (normH? normH->GetBinContent(1) : 1.0);
        
        // pu weight
        ht.fill("puwgtctr",0,plotwgts);
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
      if(chTag=="A") {
        std::vector<int> acats;
        if(isVBF)                         acats.push_back(0);
        if(isHighPt)                      acats.push_back(1);
        if(isHighPt && !isVBF && mjj>400) acats.push_back(2);
        for(auto bin : acats) ht.fill("category",bin,plotwgts);
      }
      for( auto c : chTags) {
        ht.fill(c+"_nvtx",  ev.nvtx,          plotwgts);
        ht.fill(c+"_njets", jets.size(),      plotwgts);
        ht.fill(c+"_ht",    scalarht,         plotwgts);
        ht.fill(c+"_vpt",   boson.Pt(),       plotwgts);
        ht.fill(c+"_vy",    boson.Rapidity(), plotwgts);   
        if(passJets) {
          ht.fill(c+"_vystar",     ystar, plotwgts);
          ht.fill(c+"_mjj", 	     mjj,	plotwgts);
          ht.fill(c+"_leadpt",     jets[0].Pt(),	plotwgts);
          ht.fill(c+"_subleadpt",  jets[1].Pt(),	plotwgts);
          ht.fill(c+"_centraleta", min(fabs(jets[0].Eta()),fabs(jets[1].Eta())),	plotwgts);
          ht.fill(c+"_forwardeta", max(fabs(jets[0].Eta()),fabs(jets[1].Eta())),	plotwgts);
          ht.fill(c+"_dijetpt",    (jets[0]+jets[1]).Pt(),plotwgts);
          ht.fill(c+"_detajj",     fabs(jets[0].Eta()-jets[1].Eta()),plotwgts);
          ht.fill(c+"_balance",    balance,plotwgts);

      	  ht.fill(c+"_jjetas", 	  jets[0].Eta()*jets[1].Eta(), plotwgts);
          ht.fill(c+"_dphivj0",  fabs(jets[0].DeltaPhi(boson)), plotwgts);
          ht.fill(c+"_dphivj1",  fabs(jets[1].DeltaPhi(boson)), plotwgts);
		  if(jets.size() > 2){
          	ht.fill(c+"_dphivj2",  fabs(jets[2].DeltaPhi(boson)), plotwgts);
			int nCent = 0;
		  	for(unsigned int iJet = 2; iJet < jets.size(); iJet++){	
				float dy = fabs(jets[0].Rapidity() - jets[1].Rapidity())/2;
				float sumy = (jets[0].Rapidity() + jets[1].Rapidity())/2;
				if(fabs(jets[iJet].Rapidity() - sumy) < dy){
          			ht.fill(c+"_centjy", jets[iJet].Rapidity(), plotwgts);
					nCent++;
				}
		  	}
		  	ht.fill(c+"_ncentj", nCent, plotwgts);
		  }
		  if(jets.size() > 3)
          	ht.fill(c+"_dphivj3",  fabs(jets[3].DeltaPhi(boson)), plotwgts);

		  //Study jet variables
		  for (unsigned int iJet = 0; iJet < 2; iJet++){
	  	  	ht.fill(c+"_jet_c1_00", ev.j_c1_00[jets[iJet].getJetIndex()]	  , plotwgts);
	      	ht.fill(c+"_jet_c1_02", ev.j_c1_02[jets[iJet].getJetIndex()]	  , plotwgts);
	   	    ht.fill(c+"_jet_c1_05", ev.j_c1_05[jets[iJet].getJetIndex()]	  , plotwgts);
		    ht.fill(c+"_jet_c1_10", ev.j_c1_10[jets[iJet].getJetIndex()]	  , plotwgts);
		  	ht.fill(c+"_jet_c1_20", ev.j_c1_20[jets[iJet].getJetIndex()]	  , plotwgts);

	  	  	ht.fill(c+"_jet_c2_00", ev.j_c2_00[jets[iJet].getJetIndex()]	  , plotwgts);
	      	ht.fill(c+"_jet_c2_02", ev.j_c2_02[jets[iJet].getJetIndex()]	  , plotwgts);
	  	  	ht.fill(c+"_jet_c2_05", ev.j_c2_05[jets[iJet].getJetIndex()]	  , plotwgts);
	      	ht.fill(c+"_jet_c2_10", ev.j_c2_10[jets[iJet].getJetIndex()]	  , plotwgts);
	  	  	ht.fill(c+"_jet_c2_20", ev.j_c2_20[jets[iJet].getJetIndex()]	  , plotwgts);

	  	  	ht.fill(c+"_jet_c3_00", ev.j_c3_00[jets[iJet].getJetIndex()]	  , plotwgts);
	      	ht.fill(c+"_jet_c3_02", ev.j_c3_02[jets[iJet].getJetIndex()]	  , plotwgts);
	  	  	ht.fill(c+"_jet_c3_05", ev.j_c3_05[jets[iJet].getJetIndex()]	  , plotwgts);
	      	ht.fill(c+"_jet_c3_10", ev.j_c3_10[jets[iJet].getJetIndex()]	  , plotwgts);
	  	  	ht.fill(c+"_jet_c3_20", ev.j_c3_20[jets[iJet].getJetIndex()]	  , plotwgts);

	  	  	ht.fill(c+"_jet_zg", 	ev.j_zg[jets[iJet].getJetIndex()]	  , plotwgts);
	  	  	ht.fill(c+"_jet_gaptd", ev.j_gaptd[jets[iJet].getJetIndex()]	  , plotwgts);
	  	  	ht.fill(c+"_jet_gawidth", ev.j_gawidth[jets[iJet].getJetIndex()]	  , plotwgts);
	  	  	ht.fill(c+"_jet_gathrust", ev.j_gathrust[jets[iJet].getJetIndex()]	  , plotwgts);
	  	  	ht.fill(c+"_jet_tau21", ev.j_tau21[jets[iJet].getJetIndex()]	  , plotwgts);
	  	  	ht.fill(c+"_jet_tau32", ev.j_tau32[jets[iJet].getJetIndex()]	  , plotwgts);
		  }
        }
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
