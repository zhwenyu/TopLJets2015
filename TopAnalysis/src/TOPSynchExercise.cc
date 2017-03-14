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
    allPlots["cutflow_"+fsVec[i]]  = new TH1F("cutflow_"+fsVec[i],";Cutflow;Events",9,0,9);    
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TTree *t = (TTree*)f->Get("analysis/data");
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  attachToMiniEventTree(t,ev, true);
  Int_t nentries(t->GetEntriesFast());
  cout << "...producing " << outname << " from " << nentries << " events" << endl;

  //count the original number of events processed
  Int_t ninitial=((TH1 *)f->Get("analysis/counter"))->GetBinContent(1);
  for (auto& it : allPlots) it.second->Fill(0.,ninitial);

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  SelectionTool evsel(filename,debug,triggerList);
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
            
      //decide the lepton channel
      if(ev.event==21702 || ev.event==21705 || ev.event==21927)
       	{
      	  cout << ev.run << " " << ev.lumi << " " << ev.event << " " << endl;
      	  evsel.setDebug(true);
	  cout << "MET: "<< ev.met_pt[0] << endl;
	}
      else evsel.setDebug(false);
      TString chTag = evsel.flagFinalState(ev);
      if(!evsel.passMETFilters(ev)) continue;

      //separate leptons depending on the channel being looked at
      std::vector<Particle> &leptons     = evsel.getSelLeptons();
      std::vector<Particle> &vetoLeptons = evsel.getVetoLeptons();
      std::vector<Jet>      &jets        = evsel.getJets();

      //count b-tags (CSVM by default)
      int nbtags(0);
      for (auto& jet : jets) nbtags += (abs(jet.flavor())==5);
     
      //
      if(chTag=="") continue;
      
      //veto second loose leptons in the single lepton channel
      if(chTag=="E" || chTag=="M")
	{
	  allPlots["cutflow_"+chTag]->Fill(1);

	  if(vetoLeptons.size()) continue;
	  allPlots["cutflow_"+chTag]->Fill(2);

	  if(jets.size()<1) continue;
	  allPlots["cutflow_"+chTag]->Fill(3);

	  if(nbtags>=1)
	    allPlots["cutflow_"+chTag]->Fill(4);

	  if(nbtags>=2)
	    allPlots["cutflow_"+chTag]->Fill(5);

	}
      else
	{
	  if(leptons[0].p4().Pt()<25 && leptons[1].p4().Pt()<25) continue;
	  if(leptons[0].charge()*leptons[1].charge()>0) continue;
	  allPlots["cutflow_"+chTag]->Fill(1);

	  float mll=(leptons[0].p4()+leptons[1].p4()).M();
	  if(mll<20) continue;
	  allPlots["cutflow_"+chTag]->Fill(2);
	  
	  //if(chTag=="MM")
	    //  cout << ev.run << " " << ev.lumi << " " << ev.event << " " << chTag << " " << leptons.size() << " " <<jets.size() << " " << nbtags << endl;

	  if((chTag=="EE" || chTag=="MM"))
	    {
	      if(fabs(mll-91)<15) continue;
	      allPlots["cutflow_"+chTag]->Fill(3);
	 
	      if(ev.met_pt[0]<40) continue;
	      allPlots["cutflow_"+chTag]->Fill(4);
	    }
 
	  if(jets.size()<2) continue;
	  allPlots["cutflow_"+chTag]->Fill(5);

	  if(nbtags>=1)
	    allPlots["cutflow_"+chTag]->Fill(7);

	  if(nbtags>=2)
	    allPlots["cutflow_"+chTag]->Fill(8);
	}
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

