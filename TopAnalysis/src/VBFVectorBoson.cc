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

#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"

using namespace std;

//TODOS
// pedro: jet distribution before pu id in boson+1jet
//

//
void RunVBFVectorBoson(TString filename,
                       TString outname,
                       Int_t anFlag,
                       TH1F *normH,
                       TH1F *genPU,
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

  //read photon/Z pT weights if required
  std::map<TString,TGraph *> photonPtWgts;
  std::map<TString,std::pair<double,double> > photonPtWgtCtr;
  if(anFlag>0) {
    TString wgtUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBFVectorBoson/raw/plots/ratio_plotter.root");
    gSystem->ExpandPathName(wgtUrl);
    TFile *wgtF=TFile::Open(wgtUrl);
    if(wgtF) {
      cout << "Reading photon/Z pT weights" << endl;
      TString pfix(baseName.Contains("Data13TeV_") ? "" : "_mc_MC");
      photonPtWgts["VBFA"]      = new TGraph((TH1* )wgtF->Get("VBFA_vectorbosonPt_ratio/VBFA_vectorbosonPt"+pfix));
      photonPtWgts["HighPtA"]   = new TGraph((TH1* )wgtF->Get("HighPtA_vectorbosonPt_ratio/HighPtA_vectorbosonPt"+pfix));
      photonPtWgtCtr["VBFA"]    = std::pair<double,double>(0.0,0.0);
      photonPtWgtCtr["HighPtA"] = std::pair<double,double>(0.0,0.0);
      wgtF->Close();
    } else {
      cout << "Requested to reweight photon spectrum but could not find " << wgtUrl << endl
           << "Proceeding without" << endl;
    }
  }


  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  bool isQCDEMEnriched( filename.Contains("MC13TeV_QCDEM") );

  cout << "...producing " << outname << " from " << nentries << " events" << endl;

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
  ht.addHist("puwgtctr",      new TH1F("puwgtctr",         ";Weight sums;Events",2,0,2));  
  ht.addHist("qscale",        new TH1F("qscale",           ";Q^{2} scale;Events",100,0,2000));  
  ht.addHist("nvtx",          new TH1F("nvtx",             ";Vertex multiplicity;Events",100,-0.5,99.5));
  ht.addHist("vpt", 	      new TH1F("vectorbosonPt",    ";Boson p_{T}[GeV];Events",25,0,500));
  ht.addHist("vy", 	      new TH1F("vectorbosony",     ";Boson rapidity;Events",25,0,3));
  ht.addHist("mindrl", 	      new TH1F("mindrl",           ";min #Delta R(boson,lepton);Events",25,0,6));
  ht.addHist("sihih", 	      new TH1F("sihih",            ";#sigma(i#eta,i#eta);Events",50,0,0.03));
  ht.addHist("hoe", 	      new TH1F("hoe",              ";h/e;Events",25,0,0.05));
  ht.addHist("r9", 	      new TH1F("r9",               ";r9;Events",25,0.4,1.0));
  ht.addHist("chiso", 	      new TH1F("chiso",            ";Charged isolation [GeV];Events",25,0,25));
  ht.addHist("vystar",        new TH1F("vectorbosonystar", ";y-(1/2)(y_{j1}+y_{j2});Events",25,0,5));
  ht.addHist("njets",         new TH1F("njets",            ";Jet multiplicity;Events",10,-0.5,9.5));
  ht.addHist("mjj", 	      new TH1F("mjj",              ";Dijet invariant mass [GeV];Events",40,0,4000));
  ht.addHist("leadpt",        new TH1F("leadpt",           ";Leading jet p_{T} [GeV];Events",25,0,500));
  ht.addHist("subleadpt",     new TH1F("subleadpt"   ,     ";Sub-leading jet p_{T} [GeV];Events",25,0,500));
  ht.addHist("drj1b",         new TH1F("drj1b",            ";#DeltaR(j_{1},boson);Events",25,0,8));
  ht.addHist("drj2b",         new TH1F("drj2b"   ,         ";#DeltaR(j_{2},boson);Events",25,0,8));
  ht.addHist("leadpumva",     new TH1F("leadpumva",        ";Pileup MVA;Events",25,-1,1));
  ht.addHist("subleadpumva",  new TH1F("subleadpumva"   ,  ";Pileup MVA;Events",25,-1,1));
  ht.addHist("centraleta",    new TH1F("centraleta",       ";Most central jet |#eta|;Events",25,0,5));
  ht.addHist("forwardeta",    new TH1F("forwardeta",       ";Most forward jet |#eta|;Events",25,0,5));
  ht.addHist("dijetpt",       new TH1F("dijetpt",          ";Dijet p_{T} [GeV];Events",20,0,1000));
  ht.addHist("detajj",        new TH1F("detajj" ,          ";#Delta#eta(J,J);Events",20,0,8));
  ht.addHist("dphijj",        new TH1F("dphijj" ,          ";#Delta#phi(J,J) [rad];Events",20,-3.15,3.15));
  ht.addHist("ht",            new TH1F("ht",               ";H_{T} [GeV];Events",20,0,4000));
  ht.addHist("mht",           new TH1F("mht",              ";Missing H_{T} [GeV];Events",20,0,500));
  ht.addHist("balance",       new TH1F("balance",          ";System p_{T} balance [GeV];Events",20,0,300));
  ht.addHist("relbpt",        new TH1F("relbpt",           ";#Sigma |p_{T}(j)|/Boson p_{T};Events",20,0,2));
  ht.addHist("dphibjj",       new TH1F("dphibjj",          ";#Delta#phi(JJ,boson);Events",20,-3.15,3.15));
  ht.addHist("sphericity",    new TH1F("sphericity",       ";Sphericity;Events",20,0,1.0));
  ht.addHist("aplanarity",    new TH1F("aplanarity",       ";Aplanarity;Events",20,0,1.0));
  ht.addHist("C",             new TH1F("C",                ";C;Events",20,0,1.0));
  ht.addHist("D",             new TH1F("D",                ";D;Events",20,0,1.0));
  ht.addHist("isotropy",      new TH1F("isotropy",         ";Isotropy;Events",20,0,1.0));
  ht.addHist("circularity",   new TH1F("circularity",      "Circularity;;Events",20,0,1.0));

  std::cout << "init done" << std::endl;

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, debug, triggerList,SelectionTool::VBF);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%10000==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);

      std::vector<double>plotwgts(1,1.0);
      ht.fill("qscale",ev.g_qscale,plotwgts);
      
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
      std::vector<Particle> &photons     = selector.getSelPhotons(); 
      std::vector<Particle> &leptons     = selector.getSelLeptons(); 
      std::vector<Jet>      &alljets     = selector.getJets();  
      std::vector<Jet> jets;
      for(auto j : alljets) {
        int idx=j.getJetIndex();
        int jid=ev.j_id[idx];
        bool passLoosePu((jid>>2)&0x1);
        if(!passLoosePu) continue;
        jets.push_back(j);
      }
      if(chTag!="A" && chTag!="MM") continue;

      //jet related variables and selection
      float mjj(jets.size()>=2 ?  (jets[0]+jets[1]).M() : 0.);
      float detajj(jets.size()>=2 ? fabs(jets[0].Eta()-jets[1].Eta()) : -1.);
      float dphijj(jets.size()>=2 ? jets[0].DeltaPhi(jets[1]) : -1.);
      float jjpt=(jets.size()>=2 ? (jets[0]+jets[1]).Pt() : 0.);
      float scalarht(0.);
      TLorentzVector mhtP4(0,0,0,0);
      for(auto j : jets) {
        scalarht += j.Pt();
        mhtP4 += j;
      }
      float mht(mhtP4.Pt());
      bool passJets(jets.size()>=2 && mjj>400);
      bool passVBFJetsTrigger(passJets && detajj>3.0);

      //categorize the event according to the boson kinematics
      //for the photon refine also the category according to the trigger  bit
      TLorentzVector boson(0,0,0,0);     
      bool isHighPt(false),isVBF(false),isHighPtAndVBF(false),isBosonPlusOneJet(false);
      float sihih(0),chiso(0),r9(0),hoe(0);
      if(chTag=="A") {
        
        boson += photons[0];
        sihih = ev.gamma_sieie[photons[0].originalReference()];
        chiso = ev.gamma_chargedHadronIso[photons[0].originalReference()];
        r9    = ev.gamma_r9[photons[0].originalReference()];
        hoe   = ev.gamma_hoe[photons[0].originalReference()];
        isVBF    = (selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v", ev.triggerBits) 
                    && photons[0].Pt()>75 
                    && fabs(photons[0].Eta())<1.442
                    && passVBFJetsTrigger);
        isHighPt = ( selector.hasTriggerBit("HLT_Photon200_v", ev.triggerBits) 
                     && photons[0].Pt()>200 );       
        isHighPtAndVBF = (isHighPt && isVBF);
        isBosonPlusOneJet=(isHighPt && alljets.size()==1);

        //veto prompt photons on the QCDEM enriched sample
        if( isQCDEMEnriched && ev.gamma_isPromptFinalState[ photons[0].originalReference() ] ) {
          isVBF          = false;
          isHighPt       = false;
          isHighPtAndVBF = false;
        }
          
      }
      else {
        boson   += leptons[0];
        boson   += leptons[1];
        isVBF    = boson.Pt()>75 && fabs(boson.Rapidity())<1.442 && passVBFJetsTrigger;
        isHighPt = boson.Pt()>200;
        isHighPtAndVBF = (isHighPt && isVBF);
        isBosonPlusOneJet=(isHighPt && alljets.size()==1);
      }
      if(!isVBF && !isHighPt && !isBosonPlusOneJet) continue;      

      //leptons and boson
      double mindrl(9999.);
      for(auto &l: leptons) mindrl=min(l.DeltaR(boson),mindrl);

      std::vector<TString> chTags;      
      if(isVBF)             chTags.push_back("VBF"+chTag);
      if(isHighPt)          chTags.push_back("HighPt"+chTag);
      if(isHighPtAndVBF)    chTags.push_back("HighPtVBF"+chTag);
      if(isBosonPlusOneJet) chTags.push_back("V1J"+chTag);

      //system variables and event shapes
      float ystar(0),balance(0),relbpt(0),dphibjj(0);
      if(passJets) {
        ystar=boson.Rapidity()-0.5*(jets[0].Rapidity()+jets[1].Rapidity());
        balance=(boson+jets[0]+jets[1]).Pt();
        relbpt=(jets[0].Pt()+jets[1].Pt())/boson.Pt();
        dphibjj=boson.DeltaPhi( jets[0]+jets[1] );
      }
      else if(jets.size()>0){
        balance=(boson+jets[0]).Pt();
        relbpt=jets[0].Pt()/boson.Pt();
        dphibjj=boson.DeltaPhi(jets[0]);
      }
        
      std::vector<math::XYZVector> inputVectors;
      inputVectors.push_back( math::XYZVector(boson.Px(),boson.Py(),boson.Pz()) );
      for(size_t ij=0; ij<min(size_t(2),jets.size());ij++) {
        inputVectors.push_back( math::XYZVector(jets[ij].Px(),jets[ij].Py(),jets[ij].Pz()) );
      EventShapeVariables esv(inputVectors);

      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0);
      if (!ev.isData) {

        // norm weight
        wgt  = (normH? normH->GetBinContent(1) : 1.0);
            
        // pu weight
        ht.fill("puwgtctr",0,plotwgts);
        double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
        std::vector<double>puPlotWgts(1,puWgt);
        ht.fill("puwgtctr",1,puPlotWgts);
        
        // photon trigger*selection weights
        float trigSF(1.0), selSF(1.0);
        if(chTag=="A")
          {
            trigSF *= gammaEffWR.getTriggerCorrection({},photons,{}, period).first;
            selSF  *= gammaEffWR.getOfflineCorrection(photons[0], period).first;
          }
        else
          {
            trigSF *=gammaEffWR.getTriggerCorrection(leptons,{},{}, period).first;
            selSF  *=gammaEffWR.getOfflineCorrection(leptons[0], period).first;
            selSF  *=gammaEffWR.getOfflineCorrection(leptons[1], period).first;
          }
        wgt *= puWgt*trigSF*selSF;
        
       
        // generator level weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);

        //update weight for plotter
        plotwgts[0]=wgt;
      }

      //control histograms
      for( auto c : chTags) {

        std::vector<double> cplotwgts(plotwgts);

        //photon pT weighting
        if(chTag=="A") {
          float photonPtWgt(1.0);
          if(photonPtWgts.find(c)!=photonPtWgts.end()) {
            photonPtWgt=photonPtWgts[c]->Eval(boson.Pt());
            if(photonPtWgt>0) photonPtWgt = 1./photonPtWgt;
            else              photonPtWgt = 1.0;
          }
          photonPtWgtCtr[c].first  += 1.0;
          photonPtWgtCtr[c].second += photonPtWgt;
          cplotwgts[0]*=photonPtWgt;
        } 

        ht.fill("nvtx",   ev.nvtx,          cplotwgts,c);        

        //boson histos
        ht.fill("mindrl", mindrl,           cplotwgts,c);   
        if(chTag=="A" && mindrl<0.4) continue;
        ht.fill("vpt",    boson.Pt(),       cplotwgts,c);
        ht.fill("vy",     boson.Rapidity(), cplotwgts,c);   
        ht.fill("sihih",  sihih,            cplotwgts,c);   
        ht.fill("r9",     r9,               cplotwgts,c);   
        ht.fill("hoe",    hoe,              cplotwgts,c);   
        ht.fill("chiso",  chiso,            cplotwgts,c);   

        //jet histos
        double minEta(9999),maxEta(-9999);
        for(size_t ij=0; ij<min(size_t(2),jets.size());ij++) {
          TString jtype(ij==0?"lead":"sublead");
          ht.fill(jtype+"pt",       jets[ij].Pt(),        cplotwgts,c);          
          ht.fill(jtype+"pumva",    jets[ij].getPUMVA(),  cplotwgts,c);
          ht.fill(Form("dr%db",(int)(ij+1)),   jets[ij].DeltaR(boson),  cplotwgts,c);
          minEta=min(minEta,fabs(jets[ij].Eta()));
          maxEta=max(maxEta,fabs(jets[ij].Eta()));
        }
        ht.fill("njets",        jets.size(), cplotwgts,c);
        ht.fill("ht",           scalarht,    cplotwgts,c);
        ht.fill("mht",          mht,         cplotwgts,c);
        ht.fill("centraleta",   minEta,      cplotwgts,c);
        ht.fill("forwardeta",   maxEta,      cplotwgts,c);
        ht.fill("dijetpt",      jjpt,        cplotwgts,c);
        ht.fill("detajj",       detajj,      cplotwgts,c);
        ht.fill("dphijj",       dphijj,      cplotwgts,c);
        ht.fill("mjj", 	        mjj,         cplotwgts,c);

        //visible system histos
        ht.fill("vystar",       ystar,              cplotwgts,c);        
        ht.fill("balance",      balance,            cplotwgts,c);
        ht.fill("relbpt",       relbpt,             cplotwgts,c);
        ht.fill("dphibjj",      dphibjj,            cplotwgts,c);
        ht.fill("isotropy",     esv.isotropy(),     cplotwgts,c);
        ht.fill("circularity",  esv.circularity(),  cplotwgts,c);
        ht.fill("sphericity",   esv.sphericity(1.), cplotwgts,c);
        ht.fill("aplanarity",   esv.aplanarity(1.), cplotwgts,c);
        ht.fill("C",            esv.C(1.),          cplotwgts,c);
        ht.fill("D",            esv.D(1.),          cplotwgts,c);
        }
      }
    }
  
  //close input file
  f->Close();

  //compute the scale factor needed to keep the normalization
  //due to photon pT weighting
  for(auto &wit : photonPtWgtCtr) {
    if(wit.second.second<=0) wit.second.first=1.0;
    else                     wit.second.first /= wit.second.second;
  }
  
  //save histos to file  
  fOut->cd();
  for (auto& it : ht.getPlots())  { 
    for(auto &wit : photonPtWgtCtr){
      if(!it.first.Contains(wit.first)) continue;
      cout << "Scaling " << it.first << " by "<< wit.second.first <<endl;
      it.second->Scale(wit.second.first);
      break;
    }
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    for(auto &wit : photonPtWgtCtr){
      if(!it.first.Contains(wit.first)) continue;
      it.second->Scale(wit.second.first);
      break;
    }
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}
