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
#include "TopLJets2015/TopAnalysis/interface/MttbarAnalyzer.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/TOPJetShape.h"


#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"


using namespace std;


//
void RunTTZAnalyzer(TString filename,
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
  
  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;

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
  ht.addHist("puwgtctr",     new TH1F("puwgtctr",    ";Weight sums;Events",2,0,2));
  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",55,-0.5,49.5));
  ht.addHist("njets",        new TH1F("njets",       ";Jet multiplicity;Events",15,-0.5,14.5));
  ht.addHist("nbjets",       new TH1F("nbjets",      ";b jet multiplicity;Events",10,-0.5,9.5));
  ht.addHist("mll",          new TH1F("mll",         ";M(l,l') [GeV];Events",50,91-20,91+20));
  ht.addHist("ptll",         new TH1F("ptll",        ";p_{T}(l,l') [GeV];Events",50,0,250));
  ht.addHist("dphill",       new TH1F("dphill",      ";#Delta#phi(l,l');Events",50,0,3.15));

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

      //select a Z candidate out of the leptons
      std::vector<std::pair<int,int> > zCandidates;
      for(size_t il=0; il<leptons.size(); il++)
        {
          for(size_t jl=il+1; jl<leptons.size(); jl++)
            {
              if( abs(leptons[il].id())!=abs(leptons[jl].id()) ) continue;
              if( leptons[il].charge()*leptons[jl].charge()>=0) continue;
              TLorentzVector dil=leptons[il].p4()+leptons[jl].p4();
              if( fabs(dil.M()-91)>20 ) continue;
              zCandidates.push_back( std::pair<int,int>(il,jl) );
            }
        }

      //separate jets according to flavour
      std::vector<Jet> bJets,lightJets;
      for(size_t ij=0; ij<jets.size(); ij++) 
        {
          if(jets[ij].flavor()==5) bJets.push_back(jets[ij]);
          else                     lightJets.push_back(jets[ij]);
        }

      //require at least three leptons and one Z candidate
      if(leptons.size()<3) continue;
      if(zCandidates.size()!=1) continue;
      bool passJets(lightJets.size()>=2);
      bool passBJets(bJets.size()>=2);
      if(!passJets || !passBJets) continue;
      
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
        
        // lepton trigger*selection weights
        std::vector<Particle> trigLepton(1,leptons[0]);
        EffCorrection_t trigSF = lepEffH.getTriggerCorrection(trigLepton, period);
        EffCorrection_t  selSF= lepEffH.getOfflineCorrection(leptons[0], period);

        wgt *= puWgt*trigSF.first*selSF.first;
      
        // generator level weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);

        //update weight for plotter
        plotwgts[0]=wgt;
      }

      

      //control histograms
      ht.fill("nvtx",     ev.nvtx,        plotwgts);
      ht.fill("nbjets", bJets.size(), plotwgts);
      ht.fill("njets",  jets.size(),  plotwgts);

      Particle zlepton1=leptons[ zCandidates[0].first ];
      Particle zlepton2=leptons[ zCandidates[0].second ];
      TLorentzVector z=zlepton1.p4()+zlepton2.p4();
      ht.fill("mll",     z.M(),                                  plotwgts);
      ht.fill("ptll",    z.Pt(),                                 plotwgts);
      ht.fill("dphill",  zlepton1.p4().DeltaPhi(zlepton2.p4()),  plotwgts);  
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

