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
#include "TopLJets2015/TopAnalysis/interface/BjetChargeTreeProducer.h"
#include "TopLJets2015/TopAnalysis/interface/TOPJetShape.h"


#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"


using namespace std;


//
void RunBjetChargeTreeProducer(TString filename,
                               TString outname,
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
  TTree *tree=new TTree("data","data");
  BJetSummary_t summary;
  createBJetSummaryTree(tree,summary);
  tree->SetDirectory(fOut);

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, debug);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
      
      //loop over reco jets 
      std::vector<Jet>      &jets        = selector.getJets();  
      for(size_t ij=0; ij<jets.size(); ij++)
        {
          //select only b-jets
          int origIdx=jets[ij].getJetIndex();
          if(abs(ev.j_hadflav[origIdx])!=5) continue;

          summary.pt=jets[ij].p4().Pt();
          summary.eta=jets[ij].p4().Eta();
          summary.phi=jets[ij].p4().Phi();
          summary.m=jets[ij].p4().M();
          summary.csv=jets[ij].getCSV();
          summary.nch=getMult(jets[ij]);
          summary.ptD=getPtD(jets[ij]);
          summary.ptDs=getPtDs(jets[ij]);
          summary.width=getWidth(jets[ij]);
          summary.tau21=getTau(2, 1, jets[ij]);
          summary.tau32=getTau(3, 2, jets[ij]);
          summary.tau43=getTau(4, 3, jets[ij]);
          std::vector<double> zgResult_charged = getZg(jets[ij]);
          summary.zg=zgResult_charged[0];

          summary.nch=0;
          std::vector<Particle> &pinJet=jets[ij].particles(); 
          for(size_t ipinj=0; ipinj<pinJet.size(); ipinj++)
            {
              summary.ch[summary.nch]=pinJet[ipinj].charge();
              if(summary.ch[summary.nch]==0) continue;  
              summary.chpt[summary.nch]=pinJet[ipinj].pt();
              summary.cheta[summary.nch]=pinJet[ipinj].eta();
              summary.chphi[summary.nch]=pinJet[ipinj].phi();
              summary.chm[summary.nch]=pinJet[ipinj].m();
              summary.nch++;
            }
          tree->Fill();
        }
    }
  
  //close input file
  f->Close();
  
  //save tree to file
  fOut->cd();
  tree->Write();
  fOut->Close();
}


//
void createBJetSummaryTree(TTree *t,BJetSummary_t &summary)
{
  //reco jet
  t->Branch("pt",                  &summary.pt,            "pt/F");
  t->Branch("eta",                 &summary.eta,           "eta/F");
  t->Branch("phi",                 &summary.phi,           "phi/F");
  t->Branch("m",                   &summary.m,             "m/F");
  t->Branch("csv",                 &summary.csv,           "csv/F");
  t->Branch("ptD",                 &summary.ptD,           "ptD/F");
  t->Branch("ptDs",                &summary.ptDs,         "ptDs/F");
  t->Branch("width",               &summary.width,         "width/F");
  t->Branch("tau21",               &summary.tau21,         "tau21/F");
  t->Branch("tau32",               &summary.tau32,         "tau32/F");
  t->Branch("tau43",               &summary.tau43,         "tau43/F");
  t->Branch("zg",                  &summary.zg,            "zg/F");
  t->Branch("nch",                 &summary.nch,           "nch/I");
  t->Branch("ch",                   summary.ch,            "ch[nch]/F");
  t->Branch("chpt",                 summary.chpt,          "chpt[nch]/F");
  t->Branch("cheta",                summary.cheta,         "cheta[nch]/F");
  t->Branch("chphi",                summary.chphi,         "chphi[nch]/F");
  t->Branch("chm",                  summary.chm,           "chm[nch]/F");

  //gen jet
  t->Branch("g_pt",                &summary.g_pt,            "g_pt/F");
  t->Branch("g_eta",               &summary.g_eta,           "g_eta/F");
  t->Branch("g_phi",               &summary.g_phi,           "g_phi/F");
  t->Branch("g_m",                 &summary.g_m,             "g_m/F");
  t->Branch("g_ptD",               &summary.g_ptD,           "g_ptD/F");
  t->Branch("g_ptDs",              &summary.g_ptDs,          "g_ptDs/F");
  t->Branch("g_width",             &summary.g_width,         "g_width/F");
  t->Branch("g_tau21",             &summary.g_tau21,         "g_tau21/F");
  t->Branch("g_tau32",             &summary.g_tau32,         "g_tau32/F");
  t->Branch("g_tau43",             &summary.g_tau43,         "g_tau43/F");
  t->Branch("g_zg",                &summary.g_zg,            "g_zg/F");
  t->Branch("g_bHad",              &summary.g_bHad,          "g_bHad/I");
  t->Branch("g_pId",               &summary.g_pId,           "g_pId/I");
  t->Branch("g_xb",                &summary.g_xb,            "g_xb/F");
  t->Branch("g_nch",               &summary.g_nch,           "g_nch/I");
  t->Branch("g_ch",                 summary.g_ch,            "g_ch[g_nch]/F");
  t->Branch("g_chpt",               summary.g_chpt,          "g_chpt[g_nch]/F");
  t->Branch("g_cheta",              summary.g_cheta,         "g_cheta[g_nch]/F");
  t->Branch("g_chphi",              summary.g_chphi,         "g_chphi[g_nch]/F");
  t->Branch("g_chm",                summary.g_chm,           "g_chm[g_nch]/F");
}
  
