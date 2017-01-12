#include "TopLJets2015/TopAnalysis/interface/TOPSynchExercise.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

using namespace std;


//
void RunTOPSynchExercise(TString filename, TString outname, Bool_t debug)
{

  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();

  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  std::vector<TString> fsVec = { "E", "M", "EM", "MM", "EE"};
  for(size_t i=0; i<fsVec.size(); i++)
    allPlots["cutflow_"+fsVec[i]]  = new TH1F("cutflow_"+fsVec[i],";Cutflow;Events",6,0,6);    
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev, true);
  Int_t nentries(t->GetEntriesFast());
  cout << "...producing " << outname << " from " << nentries << " events" << endl;



  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  SelectionTool evsel;
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%100==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));
      
      //decide the lepton channel
      TString chTag = evsel.flagFinalState(ev);
      if(chTag=="") continue;

      //separate leptons depending on the channel being looked at
      std::vector<Particle> &leptons     = evsel.getSelLeptons();
      std::vector<Particle> &vetoLeptons = evsel.getVetoLeptons();
      std::vector<Jet>      &jets        = evsel.getJets();

      //veto second loose leptons in the single lepton channel
      if((chTag=="E" || chTag=="M") && vetoLeptons.size()) continue;
            
      
      //for (auto& jet : jets) {
      //
      // }
    }
  
  //close input file
  f->Close();

  //save histos to file  
  fOut->cd();  
  for (auto& it : allPlots) {
    it.second->SetDirectory(fOut); 
    it.second->Write();   
  }
  fOut->Close();
}

