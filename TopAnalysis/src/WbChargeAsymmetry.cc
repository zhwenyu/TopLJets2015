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
#include "TopLJets2015/TopAnalysis/interface/WbChargeAsymmetry.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"

using namespace std;


//
void RunWbChargeAsymmetry(TString filename,
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
  WbChargeAsymmetryEvent_t tjsev;
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("wbev","wbev");
  createWbChargeAsymmetryEventTree(outT,tjsev);
  outT->SetDirectory(fOut);

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

  //lumi
  TH1F *ratevsrunH=0;
  std::map<Int_t,Float_t> lumiMap;
  if( isData )  
    {
      std::pair<std::map<Int_t,Float_t>, TH1F *> result=parseLumiInfo(era);
      lumiMap   = result.first;
      ratevsrunH = result.second;
    }
  
  //PILEUP WEIGHTING
  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if( !isData ) puWgtGr=getPileupWeightsMap(era,genPU);
  
  
  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);


  //B-TAG CALIBRATION
  BTagSFUtil* myBTagSFUtil = new BTagSFUtil();
  std::map<TString, std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> > btvsfReaders = getBTVcalibrationReadersMap(era, BTagEntry::OP_MEDIUM);
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *>    expBtagEffPy8 = readExpectedBtagEff(era);
  
  
  //JET ENERGY UNCERTAINTIES
  //std::string jecVar = "Total";
  //TString jecUncUrl(era+"/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt");
  //gSystem->ExpandPathName(jecUncUrl);
  //JetCorrectorParameters *jecParam = new JetCorrectorParameters(jecUncUrl.Data(), jecVar);
  //JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty( *jecParam );
   
  //BFRAG UNCERTAINTIES
  std::map<TString, TGraph*> bfrag = getBFragmentationWeights(era);
  std::map<TString, std::map<int, double> > semilepbr = getSemilepBRWeights(era);
  
  //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);
  std::map<TString, TH1 *> allPlots;
  std::map<TString, TH2 *> all2dPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",2,0,2);
  std::vector<TString> stageVec = { "1l", "1l1j","1l1b"};
  std::vector<TString> chTags = { "E", "M" };
  for(auto& stage : stageVec) {
    for(auto& channel : chTags) {  
      TString tag(channel+stage+"_");
      
      if(ratevsrunH) allPlots[tag+"ratevsrun"] = (TH1 *)ratevsrunH->Clone(tag+"ratevsrun");
      ht.addHist(tag+"nvtx", new TH1F(tag+"nvtx",";Vertex multiplicity;Events",55,0,55));
      ht.addHist(tag+"njets", new TH1F(tag+"njets",";Jet multiplicity;Events",15,-0.5,14.5));
      ht.addHist(tag+"nbjets", new TH1F(tag+"nbjets",";b jet multiplicity;Events",5,-0.5,4.5));
      ht.addHist(tag+"lpt", new TH1F(tag+"lpt",";Lepton p_{t} [GeV];Events",50,0,250));
      ht.addHist(tag+"leta", new TH1F(tag+"leta",";Lepton pseudo-rapidity;Events",50,-2.5,2.5));
      ht.addHist(tag+"pt", new TH1F(tag+"pt",";Jet transverse momentum [GeV];Events",50,0,250));
      ht.addHist(tag+"eta", new TH1F(tag+"eta",";Jet pseudo-rapidity;Events",50,-5,5));
      ht.addHist(tag+"met", new TH1F(tag+"met",";Missing transverse momentum [GeV];Events",50,0,250));
    }
  }
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  std::cout << "init done" << std::endl;

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, false, triggerList);
  
  for (Int_t iev=0;iev<1000+0*nentries;iev++)
    {
      t->GetEntry(iev);
      resetWbChargeAsymmetryEvent(tjsev);
      if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
      
      //assign randomly a run period
      TString period = assignRunPeriod(runPeriods,random);
      
      //////////////////
      // CORRECTIONS //
      ////////////////
      double csvm = 0.8484;
      addBTagDecisions(ev, csvm, csvm);
      
      if(!ev.isData) {
        ev = smearJetEnergies(ev);
        ev = updateBTagDecisions(ev, btvsfReaders[period],expBtagEffPy8,expBtagEffPy8,myBTagSFUtil);
      }
      
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      
      //decide the lepton channel and get selected objects
      TString chTag = selector.flagFinalState(ev);
      std::vector<Particle> &leptons     = selector.getSelLeptons(); 
      std::vector<Jet>      &jets        = selector.getJets();  
      
      //count b and W candidates
      int sel_nbjets(0);
      int seljetidx(-1);
      for(size_t ij=0; ij<jets.size(); ij++)
        {
          if (jets[ij].flavor() != 5) continue;
          ++sel_nbjets;
          if(seljetidx<0) seljetidx=ij;
        }
      
      //event selected on reco level?
      bool singleLepton         ((chTag=="E" or chTag=="M") and
                                 (selector.getVetoLeptons().size() == 0));
      bool singleLepton1Jet     (singleLepton and jets.size()>0);
      bool singleLepton1b       (singleLepton1Jet and sel_nbjets==1);
      if(!singleLepton) continue;
      
      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0);
      std::vector<double>plotwgts(1,wgt);
      allPlots["puwgtctr"]->Fill(0.,1.0);
      if (!ev.isData) {

        // norm weight
        wgt  = (normH? normH->GetBinContent(1) : 1.0);
        
        // pu weight
        double puWgt(puWgtGr[period][0]->Eval(ev.g_pu));
        allPlots["puwgtctr"]->Fill(1,puWgt);
        wgt *= puWgt;
        
        // lepton trigger*selection weights
        EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons, period);
        EffCorrection_t selSF = lepEffH.getOfflineCorrection(leptons[0], period);
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
        
        tjsev.nw = 1;
        tjsev.weight[0]=wgt;

        plotwgts[0]=wgt;
      }
      else {
        tjsev.nw=1;
        tjsev.weight[0]=1.0;
      }

      //control histograms
      TString tag(chTag);
      std::vector<TString> allTags;
      if(singleLepton)     allTags.push_back(chTag+"1l");
      if(singleLepton1Jet) allTags.push_back(chTag+"1l1j");
      if(singleLepton1b)   allTags.push_back(chTag+"1l1b");
      std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);
      for(auto &tag : allTags)
        {
          if(rIt!=lumiMap.end() && ratevsrunH) allPlots[tag+"_ratevsrun"]->Fill(std::distance(lumiMap.begin(),rIt),1./rIt->second);
          ht.fill(tag+"_nvtx", ev.nvtx, plotwgts);
          ht.fill(tag+"_njets", jets.size(), plotwgts);
          ht.fill(tag+"_nbjets", sel_nbjets, plotwgts);
          ht.fill(tag+"_lpt", leptons[0].pt(), plotwgts);
          ht.fill(tag+"_leta", leptons[0].eta(), plotwgts);
          if(jets.size())
            {
              ht.fill(tag+"_pt", jets[0].pt(), plotwgts);
              ht.fill(tag+"_eta", jets[0].eta(), plotwgts);
              ht.fill(tag+"_met", ev.met_pt[0], plotwgts);
            }
        }

      //FILL EVENT HEADER
      tjsev.run=ev.run;
      tjsev.event=ev.event;
      tjsev.lumi=ev.lumi;

      //FILL RECO TREE
      tjsev.reco_sel=singleLepton1b;
      tjsev.l_pt  = leptons[0].pt();
      tjsev.l_eta = leptons[0].eta();
      tjsev.l_phi = leptons[0].phi();
      tjsev.l_m   = leptons[0].m();
      tjsev.l_id  = leptons[0].id();
      tjsev.l_c  = leptons[0].charge();
      tjsev.met_pt=ev.met_pt[0];
      tjsev.met_phi=ev.met_phi[0];
      if(seljetidx>=0)
        {
          tjsev.j_pt      = jets[seljetidx].p4().Pt();
          tjsev.j_eta     = jets[seljetidx].p4().Eta();
          tjsev.j_phi     = jets[seljetidx].p4().Phi();
          tjsev.j_m       = jets[seljetidx].p4().M(); 
          tjsev.j_csv     = jets[seljetidx].getCSV();

          std::vector<Particle> &tks=jets[seljetidx].particles();
          tjsev.ntk=0;
          for(size_t itk=0; itk<tks.size(); itk++)
            {
              if(tks[itk].charge()==0) continue;
              tjsev.tk_pt[tjsev.ntk]=tks[itk].pt();
              tjsev.tk_eta[tjsev.ntk]=tks[itk].eta();
              tjsev.tk_phi[tjsev.ntk]=tks[itk].phi();
              tjsev.tk_c[tjsev.ntk]=tks[itk].charge();
              tjsev.tk_id[tjsev.ntk]=tks[itk].id();
              tjsev.ntk++;
            }
        }

      //GEN LEVEL TREE
      if(!isData)
        {
          //decide the lepton channel at particle level
          std::vector<Particle> genVetoLeptons = selector.getGenLeptons(ev,15.,2.5);
          std::vector<Particle> genLeptons     = selector.getGenLeptons(ev,30.,2.1);
          TString genChTag = selector.flagGenFinalState(ev, genLeptons);
          std::vector<Jet> genJets = selector.getGenJets();
     
          if(genLeptons.size())
            {
              tjsev.gl_pt=genLeptons[0].pt();
              tjsev.gl_eta=genLeptons[0].eta();
              tjsev.gl_phi=genLeptons[0].phi();
              tjsev.gl_m=genLeptons[0].m();
              tjsev.gl_id=genLeptons[0].id();
              tjsev.gl_c  = genLeptons[0].charge();
            }
   
          //count b candidates
          int sel_ngbjets(0);
          int selgjetidx(-1);
          for(size_t ij=0; ij<genJets.size(); ij++)
            {
              if(genJets[ij].flavor()!=5) continue;
              ++sel_ngbjets;
              if(selgjetidx<0) selgjetidx=ij;
            } 
         
          tjsev.ngj = genJets.size();            
          if(selgjetidx>=0)
            {
              int origgjetidx = genJets[selgjetidx].getJetIndex();
              tjsev.gj_pt     = genJets[selgjetidx].p4().Pt();
              tjsev.gj_eta    = genJets[selgjetidx].p4().Eta();
              tjsev.gj_phi    = genJets[selgjetidx].p4().Phi();
              tjsev.gj_m      = genJets[selgjetidx].p4().M();
              tjsev.gj_flavor = genJets[selgjetidx].flavor();
              tjsev.gj_bid    = ev.g_bid[origgjetidx];
              tjsev.gj_xb    = ev.g_xb[origgjetidx];
              tjsev.gj_isSemiLepBhad = ev.g_isSemiLepBhad[origgjetidx];
              std::vector<Particle> &tks=genJets[selgjetidx].particles();
              tjsev.ngtk=0;
              for(size_t itk=0; itk<tks.size(); itk++)
                {
                  if(tks[itk].charge()==0) continue;
                  tjsev.gtk_pt[tjsev.ngtk]=tks[itk].pt();
                  tjsev.gtk_eta[tjsev.ngtk]=tks[itk].eta();
                  tjsev.gtk_phi[tjsev.ngtk]=tks[itk].phi();
                  tjsev.gtk_c[tjsev.ngtk]=tks[itk].charge();
                  tjsev.gtk_id[tjsev.ngtk]=tks[itk].id();
                }
            }

          //event selected on gen level?
          bool genSingleLepton((genChTag=="E" or genChTag=="M") and
                               (genVetoLeptons.size() == 1)); 
          tjsev.gen_sel = (sel_ngbjets==1 && genSingleLepton);
        }
      
      //proceed only if event is selected on gen or reco level
      if (!tjsev.gen_sel && !tjsev.reco_sel) continue;
      
      outT->Fill();
    }
  
  //close input file
  f->Close();

  //save histos to file  
  fOut->cd();
  outT->Write();
  for (auto& it : ht.getPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : all2dPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}

//
void createWbChargeAsymmetryEventTree(TTree *t,WbChargeAsymmetryEvent_t &tjsev)
{
  //header
  t->Branch("run",      &tjsev.run,     "run/i");
  t->Branch("event",    &tjsev.event,   "event/i");
  t->Branch("lumi",     &tjsev.lumi,    "lumi/i");

  //event weights
  t->Branch("nw",       &tjsev.nw,      "nw/I");
  t->Branch("weight",    tjsev.weight,  "weight[nw]/F");

  //met
  t->Branch("met_pt",   &tjsev.met_pt,  "met_pt/F");
  t->Branch("met_phi",  &tjsev.met_phi, "met_phi/F");

  //leptons
  t->Branch("l_pt",     &tjsev.l_pt ,   "l_pt/F");
  t->Branch("l_eta",    &tjsev.l_eta ,  "l_eta/F");
  t->Branch("l_phi",    &tjsev.l_phi ,  "l_phi/F");
  t->Branch("l_m",      &tjsev.l_m ,    "l_m/F");
  t->Branch("l_id",     &tjsev.l_id ,   "l_id/I");
  t->Branch("l_c",      &tjsev.l_c ,    "l_c/I");
  t->Branch("gl_pt",    &tjsev.gl_pt ,  "gl_pt/F");
  t->Branch("gl_eta",   &tjsev.gl_eta , "gl_eta/F");
  t->Branch("gl_phi",   &tjsev.gl_phi , "gl_phi/F");
  t->Branch("gl_m",     &tjsev.gl_m ,   "gl_m/F");
  t->Branch("gl_id",    &tjsev.gl_id ,  "gl_id/I");
  t->Branch("gl_c",     &tjsev.gl_c ,   "gl_c/I");
  
  //jets
  t->Branch("j_pt",  &tjsev.j_pt ,  "j_pt/F");
  t->Branch("j_eta", &tjsev.j_eta , "j_eta/F");
  t->Branch("j_phi", &tjsev.j_phi , "j_phi/F");
  t->Branch("j_m",   &tjsev.j_m ,   "j_m/F");
  t->Branch("j_csv", &tjsev.j_csv,  "j_csv/F");
  t->Branch("gj_pt",  &tjsev.gj_pt ,  "gj_pt/F");
  t->Branch("gj_eta", &tjsev.gj_eta , "gj_eta/F");
  t->Branch("gj_phi", &tjsev.gj_phi , "gj_phi/F");
  t->Branch("gj_m",   &tjsev.gj_m ,   "gj_m/F");
  t->Branch("gj_flavor",  &tjsev.gj_flavor ,  "gj_flavor/I");
  t->Branch("gj_bid",  &tjsev.gj_bid ,  "gj_bid/I");
  t->Branch("gj_xb",   &tjsev.gj_xb ,  "gj_xb/F");
  t->Branch("gj_isSemiLepBhad",   &tjsev.gj_isSemiLepBhad ,  "gj_isSemiLepBhad/O");

  //selection flags
  t->Branch("gen_sel", &tjsev.gen_sel ,  "gen_sel/O");
  t->Branch("reco_sel", &tjsev.reco_sel ,  "reco_sel/O");
  
  //tracks associated to jets
  t->Branch("ntk",   &tjsev.ntk ,    "ntk/I");
  t->Branch("tk_c",   tjsev.tk_c ,   "tk_c[ntk]/I");
  t->Branch("tk_id",  tjsev.tk_id ,  "tk_id[ntk]/I");
  t->Branch("tk_pt",  tjsev.tk_pt ,  "tk_pt[ntk]/F");
  t->Branch("tk_eta", tjsev.tk_eta , "tk_eta[ntk]/F");
  t->Branch("tk_phi", tjsev.tk_phi , "tk_phi[ntk]/F");
  t->Branch("ngtk",   &tjsev.ngtk ,    "ngtk/I");
  t->Branch("gtk_c",   tjsev.gtk_c ,   "gtk_c[ngtk]/I");
  t->Branch("gtk_id",  tjsev.gtk_id ,  "gtk_id[ngtk]/I");
  t->Branch("gtk_pt",  tjsev.gtk_pt ,  "gtk_pt[ngtk]/F");
  t->Branch("gtk_eta", tjsev.gtk_eta , "gtk_eta[ngtk]/F");
  t->Branch("gtk_phi", tjsev.gtk_phi , "gtk_phi[ngtk]/F");
}

//
void resetWbChargeAsymmetryEvent(WbChargeAsymmetryEvent_t &tjsev)
{
  tjsev.nw=0;  
  tjsev.ntk=0;
  tjsev.ngtk=0;
  tjsev.j_pt=0;
  tjsev.l_pt=0;
  tjsev.gj_pt=0;
  tjsev.gl_pt=0;
  tjsev.gen_sel=-1;
  tjsev.reco_sel=-1;
}
