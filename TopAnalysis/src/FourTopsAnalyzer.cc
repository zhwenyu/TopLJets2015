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
#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"
#include "TopLJets2015/TopAnalysis/interface/FourTopsAnalyzer.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"

using namespace std;


//
void RunFourTopsAnalyzer(TString filename,
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
  TRandom* random = new TRandom(0); // random seed for period selection
  std::vector<RunPeriod_t> runPeriods=getRunPeriods(era);
  bool isTTbar( filename.Contains("_TTJets") or (normH and TString(normH->GetTitle()).Contains("_TTJets")));
  bool isData( filename.Contains("Data") );
  
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
  
  //PILEUP WEIGHTING
  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if( !isData ) puWgtGr=getPileupWeightsMap(era,genPU);
    
  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);

  //B-TAG CALIBRATION
  std::map<TString, std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> > btvsfReaders = getBTVcalibrationReadersMap(era, BTagEntry::OP_MEDIUM);
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *>    expBtagEffPy8 = readExpectedBtagEff(era);
  
   //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("puwgtctr", new TH1F("puwgtctr",";Weight sums;Events",2,0,2));
  ht.addHist("nvtx",     new TH1F("nvtx",";Vertex multiplicity;Events",55,-0.5,49.5));
  ht.addHist("njets",    new TH1F("njets",";Jet multiplicity;Events",15,-0.5,14.5));
  ht.addHist("nbjets",   new TH1F("nbjets",";b jet multiplicity;Events",10,-0.5,9.5));
  ht.addHist("nleptons", new TH1F("nleptons",";Lepton multiplicity;Events",6,-0.5,5.5));

  std::cout << "init done" << std::endl;

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, false, triggerList);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
      
      //assign randomly a run period
      TString period = assignRunPeriod(runPeriods,random);
      
      //////////////////
      // CORRECTIONS //
      ////////////////
      double csvm = 0.8484;
      addBTagDecisions(ev, csvm, csvm);
      if(!ev.isData) smearJetEnergies(ev);
           
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      TString chTag = selector.flagFinalState(ev);
      if(chTag=="") continue;
      std::vector<Particle> &leptons     = selector.getSelLeptons(); 
      std::vector<Jet>      &jets        = selector.getJets();  
      if(leptons.size()<2) continue;
      
      //count n b-jets
      int sel_nbjets(0);
      for(size_t ij=0; ij<jets.size(); ij++) sel_nbjets += (jets[ij].flavor()==5);

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
        double puWgt(puWgtGr[period][0]->Eval(ev.g_pu));
        std::vector<double>puPlotWgts(1,puWgt);
        ht.fill("puwgtctr",1,puPlotWgts);
        wgt *= puWgt;
        
        // lepton trigger*selection weights
        EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons, period);
        EffCorrection_t selSF(1.,0.);
        for(size_t il=0; il<leptons.size(); il++)
          { 
            EffCorrection_t  ilSelSF= lepEffH.getOfflineCorrection(leptons[il], period);
            selSF.first *= ilSelSF.first;
            selSF.second += sqrt(pow(ilSelSF.second,2)+pow(selSF.second,2));
          }
        wgt *= trigSF.first*selSF.first;
        
        //top pt weighting
        double topptsf = 1.0;
        if(isTTbar) {
          for (int igen=0; igen<ev.ngtop; igen++) {
            if(abs(ev.gtop_id[igen])!=6) continue;
            topptsf *= TMath::Exp(0.0615-0.0005*ev.gtop_pt[igen]);
          }
        }
        
        // lhe weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);

        //update weight for plotter
        plotwgts[0]=wgt;
      }

      //control histograms
      ht.fill("nvtx",     ev.nvtx,        plotwgts);
      ht.fill("njets",    jets.size(),    plotwgts);
      ht.fill("nbjets",   sel_nbjets,     plotwgts);
      ht.fill("nleptons", leptons.size(), plotwgts);
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
