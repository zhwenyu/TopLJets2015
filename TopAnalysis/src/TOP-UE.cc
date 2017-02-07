#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

//#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"
//#include "TopLJets2015/TopAnalysis/interface/TOP-16-019.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-UE.h"
#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"


#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TKey.h"

using namespace std;


//
void RunTopUE(TString filename,
	      TString outname,
	      Int_t channelSelection, 
	      Int_t chargeSelection, 
	      SelectionTool::FlavourSplitting flavourSplitting,
	      TH1F *normH, 
	      Bool_t runSysts,
	      TString era)
{
  //EVENT SELECTION WRAPPER
  SelectionTool evsel(filename,false);

  //PREPARE OUTPUT
  TopUE_t tue;
  TFile *fOut=TFile::Open(outname,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("tue","tue");
  createTopUETree(outT,tue);
  outT->SetDirectory(fOut);

  //BOOK CONTROL HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",4,0,4);
  std::vector<TString> lfsVec = { "EE", "EM", "MM" };
  for(size_t ilfs=0; ilfs<lfsVec.size(); ilfs++)   
    { 
      TString tag(lfsVec[ilfs]);
      allPlots["nvtx_"+tag]   = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events",40,0,40);
      allPlots["mll_"+tag]    = new TH1F("mll_"+tag,";Dilepton invariant mass [GeV];Events",50,0,400);
      allPlots["ptpos_"+tag]   = new TH1F("ptpos_"+tag,";Lepton transverse momentum [GeV];Events",50,20,200);
      allPlots["ptll_"+tag]    = new TH1F("ptll_"+tag,";Dilepton transverse momentum [GeV];Events",50,0,200);
      allPlots["sumpt_"+tag]   = new TH1F("sumpt_"+tag,";Transverse momentum sum [GeV];Events",50,40,300);
      allPlots["njets_"+tag]  = new TH1F("njets_"+tag,";Jet multiplicity;Events",4,2,6);
      allPlots["nbtags_"+tag] = new TH1F("nbtags_"+tag,";b-tag multiplicity;Events",5,0,5);
      allPlots["nch_"+tag]  = new TH1F("nch_"+tag,";Charged particle multiplicity;Events",200,0,200);      
      allPlots["chavgpt_"+tag]  = new TH1F("chavgpt_"+tag,";Charged particle average p_{T} [GeV];Events",100,0,20);      
      allPlots["chsumpt_"+tag]  = new TH1F("chsumpt_"+tag,";Charged particle sum p_{T} [GeV];Events",100,0,20);
    }
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }

  //check file type from name
  if(filename.Contains("Data")) runSysts=false;
  runSysts=false;
  //  bool isTTbar( filename.Contains("_TTJets") );

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      resetTopUE(tue);
      if(iev%10000==0) printf("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      //
      //RECO LEVEL analysis
      //

      //selection
      TString chTag = evsel.flagFinalState(ev);
      if(chTag=="") continue;
      std::vector<Particle> &leptons=evsel.getSelLeptons();
      TLorentzVector dil(leptons[0].p4()+leptons[1].p4());
      float mll=dil.M();
      if(mll<12) continue;
      if(leptons[0].pt()<25 && leptons[1].pt()<25)  continue;
      if(fabs(leptons[0].eta())>2.5 || fabs(leptons[1].eta())>2.5) continue;

      //divide jets
      std::vector<Jet>      &jets=evsel.getJets() ;
      std::vector<size_t> lightJetsIdx, bJetsIdx;
      for(size_t ij=0; ij<jets.size(); ij++)
	{
	  if(abs(jets[ij].flavor())==5) bJetsIdx.push_back(ij);
	  else                          lightJetsIdx.push_back(ij);
	}

      //select the charged PF candidates 
      //veto neutrals
      //veto charged within 0.01 of the selected leptons
      //veto if associated to the two leading b-jets
      std::vector<std::pair<int,TLorentzVector> > selTracks,hpTracks;
      for(int ipf = 0; ipf < ev.npf; ipf++) 
	{
	  if(ev.pf_c[ipf]==0) continue;
	  if(ev.pf_pt[ipf]<0.9) continue;
	  if(fabs(ev.pf_eta[ipf])>2.5) continue;
	
	  TLorentzVector tkP4(0,0,0,0);
	  tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);

	  bool isHP(false);
	  if(leptons.size()>0 && tkP4.DeltaR(leptons[0].p4())<0.01) isHP=true;
	  if(leptons.size()>1 && tkP4.DeltaR(leptons[1].p4())<0.01) isHP=true;
	  for(size_t ibj=0; ibj<min(bJetsIdx.size(),size_t(2)); ibj++)
	    {
	      std::vector<Particle> &pinJet=jets[ bJetsIdx[ibj] ].particles();
	      for(size_t ipinj=0; ipinj<pinJet.size(); ipinj++)
		{
		  if(pinJet[ipinj].charge()==0) continue;
		  if(pinJet[ipinj].originalReference()!=ipf) continue;
		  isHP=true;
		  break;
		}
	    }
	  
	  if(!isHP) selTracks.push_back(std::pair<int,TLorentzVector>(ipf,tkP4) );
	  else      hpTracks.push_back(std::pair<int,TLorentzVector>(ipf,tkP4) );
	}
	  
      //flag if passes selection
      tue.passSel |= (bJetsIdx.size()>=2);
      tue.nj[0]=jets.size();
      tue.nb[0]=bJetsIdx.size();

      TLorentzVector rec_tt(leptons[0].p4()+leptons[1].p4());
      if(bJetsIdx.size()>0) rec_tt += jets[ bJetsIdx[0] ].p4();
      if(bJetsIdx.size()>1) rec_tt += jets[ bJetsIdx[1] ].p4();
      rec_tt += evsel.getMET();
      tue.ptttbar[0]=rec_tt.Pt();
      tue.phittbar[0]=rec_tt.Phi();	  
    	      
      int posLepton( leptons[0].charge()>0 ? 0 : 1 );
      tue.mll[0]    = mll;
      tue.ptpos[0]  = leptons[posLepton].pt();
      tue.phipos[0] = leptons[posLepton].phi();
      tue.ptll[0]   = dil.Pt();
      tue.phill[0]  = dil.Phi();
      tue.sumpt[0]  = leptons[0].pt()+leptons[1].pt();
      tue.dphill[0] = TMath::Abs(leptons[0].p4().DeltaPhi(leptons[1].p4()));
    
      //save PF cands
      tue.n=selTracks.size();
      float nch(0.),chSumPt(0.);
      for(size_t ipf=0; ipf<selTracks.size(); ipf++)
	{
	  tue.pt[ipf]         = selTracks[ipf].second.Pt();
	  tue.eta[ipf]        = selTracks[ipf].second.Eta();
	  tue.phi[ipf]        = selTracks[ipf].second.Phi();
	  tue.isInBFlags[ipf] = 0;
	  nch++;
	  chSumPt+=tue.pt[ipf];
	}
      
      //event weight
      float wgt(1.0);
      if(!ev.isData) 
	{
	  wgt = (normH? normH->GetBinContent(1) : 1.0);
	  if(ev.g_nw>0) wgt*=ev.g_w[0];
	}
      
      
      //nominal selection control histograms
      allPlots["nvtx_"+chTag]->Fill(ev.nvtx,wgt);
      allPlots["mll_"+chTag]->Fill(tue.mll[0],wgt);
      allPlots["ptpos_"+chTag]->Fill(tue.ptpos[0],wgt);
      allPlots["ptll_"+chTag]->Fill(tue.ptll[0],wgt);
      allPlots["sumpt_"+chTag]->Fill(tue.sumpt[0],wgt);
      allPlots["nbtags_"+chTag]->Fill(tue.nb[0],wgt);
      if(tue.nb[0]>1) 
	{
	  allPlots["njets_"+chTag]->Fill(tue.nj[0],wgt);
	  allPlots["nch_"+chTag]->Fill(nch,wgt);	  
	  allPlots["chavgpt_"+chTag]->Fill(nch>0 ? chSumPt/nch : -1,wgt);
	  allPlots["chsumpt_"+chTag]->Fill(chSumPt,wgt);
	}

      //finalize ntuple
      tue.cat=0;
      if(chTag=="MM") tue.cat=13*13;
      if(chTag=="EM") tue.cat=11*13;
      if(chTag=="EE") tue.cat=11*11;
      tue.run=ev.run;
      tue.event=ev.event;
      tue.lumi=ev.lumi;
      tue.nvtx=ev.nvtx;
      tue.nw=1;
      tue.weight[0]=wgt;
   
      //all done, save it
      //outT->Fill();
    }
  
  //save histos to file  
  fOut->cd();
  outT->Write();
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}



//
void createTopUETree(TTree *t,TopUE_t &tue)
{
  //event category
  t->Branch("run",         &tue.run,         "run/I");
  t->Branch("event",       &tue.event,       "event/I");
  t->Branch("lumi",        &tue.lumi,        "lumi/I");
  t->Branch("cat",         &tue.cat,         "cat/I");
  t->Branch("gen_passSel", &tue.gen_passSel, "gen_passSel/I");
  t->Branch("passSel",     &tue.passSel,     "passSel/I");
  t->Branch("gen_nj",      &tue.gen_nj,      "gen_nj/I");
  t->Branch("gen_nb",      &tue.gen_nb,      "gen_nb/I");
  t->Branch("nj",           tue.nj,          "nj[9]/I");
  t->Branch("nb",           tue.nb,          "nb[9]/I");
  t->Branch("nvtx",        &tue.nvtx,        "nvtx/I");

  //event weights
  t->Branch("nw",     &tue.nw, "nw/I");
  t->Branch("weight",  tue.weight, "weight[nw]/F");

  //ptttbar
  t->Branch("parton_ptttbar",      &tue.parton_ptttbar,   "parton_ptttbar/F");
  t->Branch("parton_phittbar",     &tue.parton_phittbar,  "parton_phittbar/F");
  t->Branch("gen_ptttbar",         &tue.gen_ptttbar,      "gen_pttbar/F");
  t->Branch("gen_phittbar",        &tue.gen_phittbar,     "gen_phittbar/F");
  t->Branch("ptttbar",              tue.ptttbar,          "ptttbar[9]/F");
  t->Branch("phittbar",             tue.phittbar,         "phittbar[9]/F");

  //leptonic quantities
  t->Branch("ptpos",      tue.ptpos ,       "ptpos[5]/F");
  t->Branch("phipos",     tue.phipos ,      "phipos[5]/F");
  t->Branch("ptll",       tue.ptll ,        "ptll[5]/F");
  t->Branch("phill",      tue.phill ,       "phill[5]/F");
  t->Branch("mll",        tue.mll ,         "mll[5]/F");
  t->Branch("sumpt",      tue.sumpt ,       "sumpt[5]/F");
  t->Branch("dphill",     tue.dphill ,      "dphill[5]/F");
  t->Branch("gen_ptpos",  &tue.gen_ptpos ,  "gen_ptpos/F");
  t->Branch("gen_phipos", &tue.gen_phipos , "gen_phipos/F");
  t->Branch("gen_ptll",   &tue.gen_ptll ,   "gen_ptll/F");
  t->Branch("gen_phill",  &tue.gen_phill ,  "gen_phill/F");
  t->Branch("gen_mll",    &tue.gen_mll ,    "gen_mll/F");
  t->Branch("gen_sumpt",  &tue.gen_sumpt ,  "gen_sumpt/F");
  t->Branch("gen_dphill", &tue.gen_dphill , "gen_dphill/F");

  //charged particles
  t->Branch("n",          &tue.n,            "n/I");
  t->Branch("pt",          tue.pt ,          "pt[n]/F");
  t->Branch("eta",         tue.eta ,         "eta[n]/F");
  t->Branch("phi",         tue.phi ,         "phi[n]/F");
  t->Branch("isInBFlags",  tue.isInBFlags ,  "isInBFlags[n]/I");

  //gen charged particles
  t->Branch("gen_n",   &tue.gen_n,     "gen_n/I");
  t->Branch("gen_pt",   tue.gen_pt ,   "gen_pt[gen_n]/F");
  t->Branch("gen_eta",  tue.gen_eta ,  "gen_eta[gen_n]/F");
  t->Branch("gen_phi",  tue.gen_phi ,  "gen_phi[gen_n]/F");
}

//
void resetTopUE(TopUE_t &tue)
{
  //dummy event header
  tue.run=-1;  tue.lumi=0;  tue.event=0;              
  
  //reset selection flags
  tue.cat=0;   
  tue.gen_passSel=0;       
  tue.passSel=0; 

  //reset weights
  tue.nw=0;      
  
  //reset all MC truth
  tue.parton_ptttbar=0;       
  tue.parton_phittbar=0;
  tue.gen_ptttbar=0; 
  tue.gen_phittbar=0;
  tue.gen_ptpos=0;
  tue.gen_phipos=0;
  tue.gen_ptll=0;
  tue.gen_mll=0;
  tue.gen_sumpt=0;
  tue.gen_dphill=0;
  tue.gen_nj=0;             
  tue.gen_nb=0;
  
  //reset particle counters
  tue.n=0; 
  tue.gen_n=0; 
}
