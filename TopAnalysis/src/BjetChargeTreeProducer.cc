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
#include "TopQuarkAnalysis/BFragmentationAnalyzer/interface/BFragmentationAnalyzerUtils.h"

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
      selector.flagFinalState(ev);
      if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);

      //loop over reco jets
      std::vector<Jet>      &jets        = selector.getJets();  
      for(size_t ij=0; ij<jets.size(); ij++)
        {

          //require an associated gen jet
          int genJetIdx=ev.j_g[ij];
          if(genJetIdx<0) continue;

          //require that there was a B-hadron matched to this jet
          summary.g_bId=ev.g_bid[genJetIdx];
          if( !IS_BHADRON_PDGID(summary.g_bId) ) continue;

          //save kinematics for this jet
          summary.pt        = jets[ij].p4().Pt();
          summary.eta       = jets[ij].p4().Eta();
          summary.phi       = jets[ij].p4().Phi();
          summary.m         = jets[ij].p4().M();

          //save kinematics of the generator-level jet
          summary.g_pt=ev.g_pt[genJetIdx];
          summary.g_eta=ev.g_eta[genJetIdx];
          summary.g_phi=ev.g_phi[genJetIdx];
          summary.g_m=ev.g_m[genJetIdx];
          summary.g_bId=ev.g_bid[genJetIdx];
          summary.g_pId=ev.g_id[genJetIdx];
          summary.g_xb=ev.g_xb[genJetIdx];
          
          //secondary vertex information (if available)
          summary.vtxmass   = ev.j_vtxmass[ jets[ij].getJetIndex() ];
          summary.vtxchi2   = ev.j_vtxchi2[ jets[ij].getJetIndex() ];
          summary.vtxL3d    = ev.j_vtx3DVal[ jets[ij].getJetIndex() ];
          summary.vtxL3dSig = ev.j_vtx3DSig[ jets[ij].getJetIndex() ];
          
          //charge estimators
          float alpha=1.1;
          std::vector<Particle> &pinJet=jets[ij].particles(); 
          summary.nch=0;
          summary.vtxnch=0;
          summary.nmu=0;
          summary.much=0;
          summary.vtxnmu=0;
          summary.vtxmuch=0;
          float avgch(0.), avgch_wgts(0.);
          float avgchVtx(0.), avgchVtx_wgts(0.);
          for(size_t ipinj=0; ipinj<pinJet.size(); ipinj++)
            {

              //require particle to have charge
              if(pinJet[ipinj].charge()==0) continue;

              summary.nch++;
              float chwgt = pow(pinJet[ipinj].pt(),alpha);
              avgch      += pinJet[ipinj].charge()*chwgt;
              avgch_wgts += chwgt;

              //check it it's a muon
              if( abs(pinJet[ipinj].id())==13 )
                {
                  summary.nmu++;
                  summary.much+=pinJet[ipinj].charge();
                }

              //check if this track belongs to a secondary vertex
              if(pinJet[ipinj].qualityFlags()==0) continue;
              summary.vtxnch++;
              avgchVtx      += pinJet[ipinj].charge()*chwgt;
              avgchVtx_wgts += chwgt;

              //check if it's a muon in the secondary vertex
              if( abs(pinJet[ipinj].id())==13 )
                {
                  summary.vtxnmu++;
                  summary.vtxmuch+=pinJet[ipinj].charge();
                }
            }
          //finalize the charge estimators
          summary.ch    = avgch_wgts==0? 0 : avgch/avgch_wgts;
          summary.vtxch = avgchVtx_wgts==0 ? 0 : avgchVtx/avgchVtx_wgts;

          //fill the tree
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
  t->Branch("vtxmass",             &summary.vtxmass,       "vtxmass/F");
  t->Branch("vtxchi2",             &summary.vtxchi2,       "vtxchi2/F");
  t->Branch("vtxL3d",              &summary.vtxL3d,        "vtxL3d/F");
  t->Branch("vtxL3dSig",           &summary.vtxL3dSig,     "vtxL3dSig/F");
  t->Branch("nch",                 &summary.nch,           "nch/I");
  t->Branch("ch",                  &summary.ch,            "ch/F");
  t->Branch("vtxnch",              &summary.vtxnch,        "vtxnch/I");
  t->Branch("vtxch",               &summary.vtxch,         "vtxch/F");
  t->Branch("nmu",                 &summary.nmu,           "nmu/I");
  t->Branch("much",                &summary.much,          "much/F");
  t->Branch("vtxnmu",              &summary.vtxnmu,        "vtxnmu/I");
  t->Branch("vtxmuch",             &summary.vtxmuch,       "vtxmuch/F");
  
  //gen jet
  t->Branch("g_pt",                &summary.g_pt,            "g_pt/F");
  t->Branch("g_eta",               &summary.g_eta,           "g_eta/F");
  t->Branch("g_phi",               &summary.g_phi,           "g_phi/F");
  t->Branch("g_m",                 &summary.g_m,             "g_m/F");
  t->Branch("g_xb",                &summary.g_xb,            "g_xb/F");
  t->Branch("g_bId",               &summary.g_bId,           "g_bId/I");
  t->Branch("g_pId",               &summary.g_pId,           "g_pId/I");
}
  
