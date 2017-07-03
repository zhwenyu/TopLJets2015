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
          float avgchPow[]={0.5,0.8,1.0,1.5};
          std::vector<float> avgch(4,0.),avgch2(4,0.),avgch3(4,0.),avgchwgts(4,0.); 
          for(size_t ipinj=0; ipinj<pinJet.size(); ipinj++)
            {
              summary.ch[summary.nch]=pinJet[ipinj].charge();
              if(summary.ch[summary.nch]==0) continue;  
              summary.chpt[summary.nch]=pinJet[ipinj].pt();
              summary.cheta[summary.nch]=pinJet[ipinj].eta();
              summary.chphi[summary.nch]=pinJet[ipinj].phi();
              summary.chm[summary.nch]=pinJet[ipinj].m();


              for(size_t ipow=0; ipow<4; ipow++)
                {
                  float chwgt=pow(summary.chpt[summary.nch],avgchPow[ipow]);
                  avgch[ipow]  += summary.ch[summary.nch]*chwgt;
                  avgch2[ipow] += pow(summary.ch[summary.nch],2)*chwgt;
                  avgch3[ipow] += pow(summary.ch[summary.nch],3)*chwgt;
                  avgchwgts[ipow]    += chwgt;
                }

              summary.ch_05 = avgch[0]/avgchwgts[0];
              summary.ch_08 = avgch[1]/avgchwgts[1];
              summary.ch_1  = avgch[2]/avgchwgts[2];
              summary.ch_15 = avgch[3]/avgchwgts[3];
              summary.ch2_05 = avgch2[0]/avgchwgts[0];
              summary.ch2_08 = avgch2[1]/avgchwgts[1];
              summary.ch2_1  = avgch2[2]/avgchwgts[2];
              summary.ch2_15 = avgch2[3]/avgchwgts[3];
              summary.ch3_05 = avgch3[0]/avgchwgts[0];
              summary.ch3_08 = avgch3[1]/avgchwgts[1];
              summary.ch3_1  = avgch3[2]/avgchwgts[2];
              summary.ch3_15 = avgch3[3]/avgchwgts[3];

              //gen jet
              int genJetIdx=ev.j_g[ij];
              if(genJetIdx>=0)
                {
                  summary.g_pt=ev.g_pt[genJetIdx];
                  summary.g_eta=ev.g_eta[genJetIdx];
                  summary.g_phi=ev.g_phi[genJetIdx];
                  summary.g_m=ev.g_m[genJetIdx];
                  summary.g_bId=ev.g_bid[genJetIdx];
                  summary.g_pId=ev.g_id[genJetIdx];
                  summary.g_xb=ev.g_xb[genJetIdx];
                }
              else
                {
                  summary.g_pt=0;
                  summary.g_eta=0;
                  summary.g_phi=0;
                  summary.g_m=0;
                  summary.g_bId=0;
                  summary.g_pId=0;
                  summary.g_xb=0;
                }

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
  t->Branch("ch_05",               &summary.ch_05,         "ch_05/F");
  t->Branch("ch_08",               &summary.ch_08,         "ch_08/F");
  t->Branch("ch_1",                &summary.ch_1,          "ch_1/F");
  t->Branch("ch_15",               &summary.ch_15,         "ch_15/F");
  t->Branch("ch2_05",              &summary.ch2_05,        "ch2_05/F");
  t->Branch("ch2_08",              &summary.ch2_08,        "ch2_08/F");
  t->Branch("ch2_1",               &summary.ch2_1,         "ch2_1/F");
  t->Branch("ch3_15",              &summary.ch2_15,        "ch2_15/F");
  t->Branch("ch3_05",              &summary.ch3_05,        "ch3_05/F");
  t->Branch("ch3_08",              &summary.ch3_08,        "ch3_08/F");
  t->Branch("ch3_1",               &summary.ch3_1,         "ch3_1/F");
  t->Branch("ch3_15",              &summary.ch3_15,        "ch3_15/F");

  //gen jet
  t->Branch("g_pt",                &summary.g_pt,            "g_pt/F");
  t->Branch("g_eta",               &summary.g_eta,           "g_eta/F");
  t->Branch("g_phi",               &summary.g_phi,           "g_phi/F");
  t->Branch("g_m",                 &summary.g_m,             "g_m/F");
  t->Branch("g_xb",                &summary.g_xb,            "g_xb/F");
  t->Branch("g_bId",               &summary.g_bId,           "g_bId/I");
  t->Branch("g_pId",               &summary.g_pId,           "g_pId/I");
}
  
