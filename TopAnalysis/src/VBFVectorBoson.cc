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
  ht.addHist("puwgtctr",     new TH1F("puwgtctr",      ";Weight sums;Events",2,0,2));
  ht.addHist("nvtx",         new TH1F("nvtx",          ";Vertex multiplicity;Events",100,-0.5,101.5));
  ht.addHist("njets",        new TH1F("njets",         ";Jet multiplicity;Events",15,-0.5,14.5));
  ht.addHist("mjj", 	     new TH1F("dijet_mass",    ";M_{jj} [GeV];Events",80,400,2000));
  ht.addHist("vpt", 	     new TH1F("vectorbosonPt", "p_{T}(V) [GeV];Events",25,50,500));
  ht.addHist("leadjetpt",    new TH1F("leadjetPt",     "p_{T}(lead. J) [GeV];Events",25,50,500));
  ht.addHist("subleadjetpt", new TH1F("subleadjetPt",  "p_{T}(sublead. J) [GeV];Events",25,50,500));
  ht.addHist("dijetpt",      new TH1F("dijetPt",       "p_{T}(JJ) [GeV];Events",25,50,500));
  ht.addHist("detajj", 	     new TH1F("DEta_jj",       "#Delta#eta(J,J);Events",20,-5,5));

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
      if(chTag=="") continue;
      std::vector<Particle> &photons     = selector.getSelPhotons(); 
      std::vector<Jet>      &jets        = selector.getJets();  


      //require one good photon
	  std::cout<<"Photon size: "<<photons.size()<<endl;
	  std::cout<<"Jet size: "<<jets.size()<<endl;	
      if(photons.size()!=1) continue;
      bool passJets(jets.size()>=2);

      
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
        EffCorrection_t trigSF = gammaEffWR.getTriggerCorrection({},photons,{}, period);
        EffCorrection_t  selSF= gammaEffWR.getOfflineCorrection(photons[0], period);

        wgt *= puWgt*trigSF.first*selSF.first;
        
        
        // generator level weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);

        //update weight for plotter
        plotwgts[0]=wgt;
      }



      //control histograms
      ht.fill("nvtx",     ev.nvtx,        plotwgts);
	  ht.fill("njets",    jets.size(),    plotwgts);
      if(passJets){
  		ht.fill("mjj", 	(jets[0]+jets[1]).M(),	plotwgts);
  		ht.fill("vpt", 	photons[0].Pt(),		plotwgts);
  		ht.fill("leadjetpt", jets[0].Pt(),		plotwgts);
  		ht.fill("subleadjetpt", jets[1].Pt(),	plotwgts);
  		ht.fill("dijetpt",(jets[0]+jets[1]).Pt(),plotwgts);
  		ht.fill("detajj", jets[0].Eta()-jets[1].Eta(),plotwgts);
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
