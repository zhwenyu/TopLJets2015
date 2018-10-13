#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/ExclusiveTop.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include <string>

#include "TMath.h"

#include "TopLJets2015/CTPPSAnalysisTools/interface/LHCConditionsFactory.h"
#include "TopLJets2015/CTPPSAnalysisTools/interface/AlignmentsFactory.h"
#include "TopLJets2015/CTPPSAnalysisTools/interface/XiReconstructor.h"

using namespace std;

#define ADDVAR(x,name,t,tree) tree->Branch(name,x,TString(name)+TString(t))

//TODO
//check PPS code is the latest from Laurent
//PPS json
//lumi, puweighting, genweighting
//launch first ntuplization on condor

//
void RunExclusiveTop(TString filename,
                     TString outname,
                     Int_t channelSelection, 
                     Int_t chargeSelection, 
                     TH1F *normH, 
                     TH1F *genPU, 
                     TString era,
                     Bool_t debug)
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  const char* CMSSW_BASE = getenv("CMSSW_BASE");
  // CTPPS reconstruction part
  ctpps::XiReconstructor proton_reco;
  proton_reco.feedDispersions(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/dispersions.txt", CMSSW_BASE));

  ctpps::AlignmentsFactory ctpps_aligns;
  ctpps_aligns.feedAlignments(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/alignments_21aug2018.txt", CMSSW_BASE));

  ctpps::LHCConditionsFactory lhc_conds;
  lhc_conds.feedConditions(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/xangle_tillTS2.csv", CMSSW_BASE));
  lhc_conds.feedConditions(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/xangle_afterTS2.csv", CMSSW_BASE));

  bool isTTbar( filename.Contains("_TTJets") or (normH and TString(normH->GetTitle()).Contains("_TTJets")));
  
  MiniEvent_t ev;
  Int_t evcat(0);


  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("tree","tree");
  outT->Branch("run",&ev.run,"run/i");
  outT->Branch("event",&ev.event,"event/l");
  outT->Branch("lumi",&ev.lumi,"lumi/i");
  outT->Branch("nvtx",&ev.nvtx,"nvtx/I");
  outT->Branch("metfilters",&ev.met_filterBits,"metfilters/I");  
  ADDVAR(&ev.rho,"rho","F",outT);
  ADDVAR(&ev.pf_closestDZnonAssoc,"closestDZnonAssoc","F",outT);
  ADDVAR(&ev.pf_ch_wgtSum,"nch","F",outT);
  ADDVAR(&ev.met_pt[1],"met_pt","F",outT);
  ADDVAR(&ev.met_phi[1],"met_phi","F",outT);
  ADDVAR(&ev.met_sig[1],"met_sig","F",outT);
  TString fvars[]={"evwgt", "evcat",
                   "l1pt", "l1eta", "l1phi", "ml1", "l1id", "mt1",
                   "l2pt", "l2eta", "l2phi", "ml2", "l2id", "mt2",
                   "llpt", "lleta", "llphi", "mll", "llht", "llacopl", "llcosthetaCS", "llphistar", "llMR", "llR", "mtll", "llcsip", "llcsim",
                   "xpt",  "xeta",  "xphi",  "mx",  "xid", "xht", "mtx",
                   "fpt",  "feta",  "fphi",  "mf",  "fht", "facopl", "fcosthetaCS", "fphistar", "fMR","fR", "fcsip", "fcsim",
                   "nb", "nj", "nl","ng","nch", "ht", "htb", "htj", "closestkdz",
                   "puppirecoil_pt","puppirecoil_phi", "puppirecoil_spher",
  };
  std::map<TString,Float_t> outVars;
  for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++){
    outVars[fvars[i]]=0.;
    ADDVAR(&(outVars[fvars[i]]),fvars[i],"F",outT);
  }
  int nRPtk(0),RPid[200];
  float RPcsi[200],RPcsiUnc[200],RPxang[200];
  if(filename.Contains("Data13TeV")){
    outT->Branch("nRPtk",&nRPtk,"nRPtk/i");
    outT->Branch("RPid",RPid,"RPid[nRPtk]/i");
    outT->Branch("RPcsi",RPcsi,"RPcsi[nRPtk]/F");
    outT->Branch("RPcsiUnc",RPcsiUnc,"RPcsiUnc[nRPtk]/F");
    outT->Branch("RPxang",RPxang,"RPxang[nRPtk]/F");
  }
  outT->SetDirectory(fOut);

  //READ TREE FROM FILE
  TFile *f = TFile::Open(filename);  
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
  //LUMINOSITY+PILEUP
  LumiTools lumi(era,genPU);
    
  //LEPTON EFFICIENCIES
  EfficiencyScaleFactorsWrapper lepEffH(filename.Contains("Data13TeV"),era);

  //B-TAG CALIBRATION
  BTagSFUtil btvSF(era,"DeepCSV",BTagEntry::OperatingPoint::OP_MEDIUM,"",0);
  
   //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("puwgtctr",     new TH1F("puwgtctr",    ";Weight sums;Events",2,0,2));
  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",50,-0.5,49.5));
  ht.addHist("nlep",         new TH1F("nlep",        ";Lepton multipliciy;Events",15,2,5));
  ht.addHist("nljets",       new TH1F("nljets",      ";light jet multiplicity;Events",5,0,5)); 
  ht.addHist("nbjets",       new TH1F("nbjets",      ";b jet multiplicity;Events",5,0,5));
  ht.addHist("lmpt",         new TH1F("lmpt",        ";Lepton 1 transverse momentum [GeV]",50,20,200));
  ht.addHist("lmeta",        new TH1F("lmeta",       ";Lepton 1 pseudo-rapidity",10,0,2.5));
  ht.addHist("lppt",         new TH1F("lppt",        ";Lepton 2 transverse momentum [GeV]",50,20,200));
  ht.addHist("lpeta",        new TH1F("lpeta",       ";Lepton 2 pseudo-rapidity",10,0,2.5));
  ht.addHist("mll",          new TH1F("mll",         ";Dilepton invariant mass [GeV]",50,20,200));
  ht.addHist("ptll",         new TH1F("ptll",        ";Dilepton transverse momentum [GeV]",50,0,200));  
  ht.addHist("phistar",      new TH1F("phistar",     ";log_{10}(dilepton #phi^{*})",20,-3,1));
  ht.addHist("costhetaCS",   new TH1F("costhetaCS",  ";Dilepton cos#theta^{*}_{CS}",20,-1,1));
  ht.addHist("acopl",        new TH1F("acopl",       ";Acoplanarity",20,0,0.25));

  std::cout << "init done" << std::endl;
  if (debug){std::cout<<"\n DEBUG MODE"<<std::endl;}

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, false, triggerList);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
	
      //assign randomly a run period
      TString period = lumi.assignRunPeriod();

      //////////////////
      // CORRECTIONS //
      ////////////////      
      btvSF.addBTagDecisions(ev);
      btvSF.updateBTagDecisions(ev);      
           
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      //trigger
      bool hasETrigger(selector.hasTriggerBit("HLT_Ele32_eta2p1_WPTight_Gsf_v", ev.triggerBits) ||
                       selector.hasTriggerBit("HLT_Ele35_eta2p1_WPTight_Gsf_v", ev.triggerBits) ||
                       selector.hasTriggerBit("HLT_Ele38_eta2p1_WPTight_Gsf_v", ev.triggerBits) ||
                       selector.hasTriggerBit("HLT_Ele40_eta2p1_WPTight_Gsf_v", ev.triggerBits) );
      bool hasMTrigger(selector.hasTriggerBit("HLT_IsoMu24_2p1_v", ev.triggerBits) ||
                       selector.hasTriggerBit("HLT_IsoMu27_v",     ev.triggerBits) );
      bool hasMMTrigger(selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                  ev.triggerBits) ||
                        selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",          ev.triggerBits) ||
                        selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",        ev.triggerBits) );
      bool hasEETrigger(selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",             ev.triggerBits) ||
                        selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          ev.triggerBits) );
      bool hasEMTrigger(selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev.triggerBits) ||
                        selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits) ||
                        selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",    ev.triggerBits) ||
                        selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits) ||
                        selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     ev.triggerBits) ||
                        selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev.triggerBits) );

      //trigger specific to data
      if (ev.isData) {

        //use only these unprescaled triggers for these eras
        if(filename.Contains("2017E") || filename.Contains("2017F")){
          hasMTrigger=selector.hasTriggerBit("HLT_IsoMu27_v",ev.triggerBits);
        }
        if(!(filename.Contains("2017A") || filename.Contains("2017B"))){
          hasMMTrigger=(selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",   ev.triggerBits) ||
                        selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev.triggerBits) );
        }

        //remove trigger overlaps (use single lepton PDs only for single lepton triggers)
        if(selector.isDoubleEGPD()       && !hasEETrigger) continue;
        if(selector.isDoubleMuonPD()     && !hasMMTrigger) continue;
        if(selector.isMuonEGPD()         && !hasEMTrigger) continue;
        if(selector.isSingleElectronPD() && (!hasETrigger || hasEETrigger || hasEMTrigger) ) continue;
        if(selector.isSingleMuonPD()     && (!hasMTrigger || hasMMTrigger || hasEMTrigger) ) continue;
      }

      std::vector<Particle> flaggedLeptons = selector.flaggedLeptons(ev);
      std::vector<Particle> leptons        = selector.selLeptons(flaggedLeptons,SelectionTool::TIGHT);
      if(leptons.size()<2) continue;
      if(leptons[0].Pt()<30 || fabs(leptons[0].Eta())>2.1) continue;

      //Z finder (leading pT leptons)
      int l1idx(0),l2idx(1);
      for(size_t il=0; il<leptons.size(); il++)
        for(size_t jl=1; jl<leptons.size(); jl++) {
          if(abs(leptons[il].id())!=abs(leptons[jl].id())) continue;
          float mll((leptons[il]+leptons[jl]).M());
          if( fabs(mll-90)>10 ) continue;
          l1idx=il;
          l2idx=jl;
          break;
        }
      bool isSF( leptons[l1idx].id()==leptons[l2idx].id() );
      bool isSS( leptons[l1idx].charge()*leptons[l2idx].charge() > 0 );

      TLorentzVector lm(leptons[l1idx].charge()>0 ? leptons[l1idx] : leptons[l1idx]);
      TLorentzVector lp(leptons[l1idx].charge()>0 ? leptons[l2idx] : leptons[l2idx]);
      if(isSS)  { lm=leptons[l1idx]; lp=leptons[l2idx]; }
      TLorentzVector dil(lm+lp);
      float mll(dil.M());
      bool isZ( isSF && !isSS && fabs(mll-91)<10);
 
      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt[1],0,ev.met_phi[1],0.);

      //photons
      std::vector<Particle> flaggedPhotons = selector.flaggedPhotons(ev);
      std::vector<Particle> photons        = selector.selPhotons(flaggedPhotons,SelectionTool::TIGHT,leptons,30.,2.5);
      
      //jets (require PU jet id)
      std::vector<Jet> allJets = selector.getGoodJets(ev,25.,2.5,leptons,photons);
      int njets(0);
      std::vector<Jet> bJets,lightJets,jets;
      float scalarht(0.),scalarhtb(0.),scalarhtj(0.);
      for(size_t ij=0; ij<allJets.size(); ij++) 
        {
          int idx=allJets[ij].getJetIndex();
          int jid=ev.j_id[idx];
          bool passLoosePu((jid>>2)&0x1);
          if(!passLoosePu) continue;
          if(allJets[ij].flavor()==5) { bJets.push_back(allJets[ij]);     scalarhtb+=allJets[ij].pt(); }
          else                        { lightJets.push_back(allJets[ij]); scalarhtj+= allJets[ij].pt(); }
          njets++;
          jets.push_back(allJets[ij]);
          scalarht += jets[ij].pt();
        }


      ///dilepton kinematics
      Float_t llht(lm.Pt()+lp.Pt());
      Float_t llacopl = computeAcoplanarity(lm,lp);
      Float_t llphistar=computePhiStar(lm,lp);
      Float_t llcosthetaCS=computeCosThetaStar(lm,lp);
      Float_t llMR(computeMR(lm,lp));
      Float_t llR(computeRsq(lm,lp,met));
      Float_t mtm(computeMT(lm,met));
      Float_t mtp(computeMT(lp,met));

      //baseline categories and additional stuff produced with the dilepton
      std::vector<TString> cats(1,"inc");

      TString dilCat(isSS ? "ss" : "os" );
      dilCat+=isSF ? "sf" : "of";
      if(isZ) dilCat += "z";
      cats.push_back(dilCat);       
      evcat=(isSS | isSF<<1 | isZ<<2);

      TLorentzVector X(met);
      Int_t xid(0);
      Float_t xht(0);
      evcat |= (1<<7);
      cats.push_back(dilCat+"inv");
      if(leptons.size()>2){
        evcat |= (1<<3);
        cats[2]=dilCat+"3l"; 
        int l3idx(1);
        if(l1idx==0) {
          if(l2idx==1) l3idx=2;
        }
        else if(l1idx==1) {
          l3idx=0;
        }
        X=leptons[l3idx];
        xid=leptons[l3idx].id();
        xht=X.Pt();
      }
      else if(photons.size()>0) {
        evcat |= (1<<4);
        cats[2]=dilCat+"a";
        X=photons[0].p4();
        xid=22;
        xht=X.Pt();
      }
      else if (jets.size()>2) {
        evcat |= (1<<5);
        cats[2]=dilCat+"jj";
        X=jets[0]+jets[1];
        xid=2121;
        xht=jets[0].Pt()+jets[1].Pt();
        if(bJets.size()>1) {
          evcat |= (1<<6);
          cats[2]=dilCat+"bb";
          X=bJets[0]+bJets[1];
          xid=55;
          xht=bJets[0].Pt()+bJets[1].Pt();
        }
      }

      //final state F=ll+X
      TLorentzVector F(dil+X);
      Float_t fht(llht+xht);
      Float_t facopl(computeAcoplanarity(dil,X));
      Float_t fphiStar(computePhiStar(dil,X));
      Float_t fcosthetaStarCS(computeCosThetaStar(dil,X));
      Float_t fMR(computeMR(dil,X));
      Float_t fR(computeRsq(dil,X,met));
      Float_t mtll(computeMT(dil,met));
      Float_t mtx(computeMT(X,met));

      //recoil and UE
      int nch(int(ev.pf_ch_wgtSum));
      float closestTrackDZ(ev.pf_closestDZnonAssoc);
      TVector2 puppiRecoil(ev.pf_puppi_px,ev.pf_puppi_py);
      float puppiRecoilHt(ev.pf_puppi_ht);
      TVector2 h2( xid!=0 ? TVector2(F.Px(),F.Py()) : TVector2(dil.Px(),dil.Py()) );
      puppiRecoil   -= h2;
      puppiRecoilHt -= (xid!=0 ? fht : llht);


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
        double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
        std::vector<double>puPlotWgts(1,puWgt);
        ht.fill("puwgtctr",1,puPlotWgts);
        
        // lepton trigger*selection weights
        EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons,{},{},period);
        EffCorrection_t  selSF = lepEffH.getOfflineCorrection(leptons[0], period);

        wgt *= puWgt*trigSF.first*selSF.first;
        
        //top pt weighting
        double topptsf = 1.0;
        if(isTTbar) {
          for (int igen=0; igen<ev.ngtop; igen++) {
            if(abs(ev.gtop_id[igen])!=6) continue;
            topptsf *= TMath::Exp(0.0615-0.0005*ev.gtop_pt[igen]);
          }
        }
        
        // generator level weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);

        //update weight for plotter
        plotwgts[0]=wgt;
      }

      //control histograms
      ht.fill("nvtx",       ev.nvtx,         plotwgts, cats);
      
      //dilepton system
      ht.fill("nlep",       leptons.size(),  plotwgts, cats);
      ht.fill("lmpt",       lm.Pt(),         plotwgts, cats);
      ht.fill("lmeta",      fabs(lm.Eta()),  plotwgts, cats);
      ht.fill("lppt",       lp.Pt(),         plotwgts, cats);
      ht.fill("lmeta",      fabs(lp.Eta()),  plotwgts, cats);
      ht.fill("mll",        dil.M(),         plotwgts, cats);
      ht.fill("ptll",       dil.Pt(),        plotwgts, cats);
      ht.fill("phistar",    TMath::Log10(llphistar),  plotwgts, cats);
      ht.fill("costhetaCS", llcosthetaCS,      plotwgts, cats);
      ht.fill("acopl",      llacopl,           plotwgts, cats);

      //bjets
      ht.fill("nbjets",     bJets.size(),    plotwgts, cats);
      ht.fill("scalarhtb",     scalarhtb,    plotwgts, cats);
      for(auto &b:bJets) {
        ht.fill("bpt",     b.Pt(),           plotwgts, cats);
        ht.fill("beta",    fabs(b.Eta()),    plotwgts, cats);
      }

      //light jets
      ht.fill("nljets",     lightJets.size(),    plotwgts, cats);
      ht.fill("scalarhtj",  scalarhtj,    plotwgts, cats);
      for(auto &l:lightJets) {
        ht.fill("jpt",     l.Pt(),           plotwgts, cats);
        ht.fill("jeta",    fabs(l.Eta()),    plotwgts, cats);
      }

      //photons
      ht.fill("npho",       photons.size(),  plotwgts, cats);
      for(auto &a:photons) {
        ht.fill("apt",     a.Pt(),           plotwgts, cats);
        ht.fill("aeta",    fabs(a.Eta()),    plotwgts, cats);
      }

      ht.fill("met",           met.Pt(),                         plotwgts, cats);
      ht.fill("nch",           nch,                              plotwgts, cats);

      //fill tree
      outVars["evcat"]=float(evcat);
      outVars["evwgt"]=plotwgts[0];

      outVars["l1pt"]=lm.Pt();
      outVars["l1eta"]=lm.Eta();
      outVars["l1phi"]=lm.Phi(); 
      outVars["ml1"]=lm.M();
      outVars["l1id"]=leptons[l1idx].id();
      outVars["mt1"]=mtm;

      outVars["l2pt"]=lp.Pt();
      outVars["l2eta"]=lp.Eta();
      outVars["l2phi"]=lp.Phi();
      outVars["ml2"]=lp.M();
      outVars["l2id"]=leptons[l2idx].id();
      outVars["mt2"]=mtp;

      outVars["llpt"]=dil.Pt();
      outVars["lleta"]=dil.Eta();
      outVars["llphi"]=dil.Phi();
      outVars["mll"]=dil.M();
      outVars["llacopl"]=llacopl;
      outVars["llcosthetaCS"]=llcosthetaCS;
      outVars["llphistar"]=llphistar;
      outVars["llMR"]=llMR;
      outVars["llR"]=llR;
      outVars["mtll"]=mtll;
      outVars["llht"]=llht;      

      std::pair<Float_t,Float_t> llcsi=calcCsi(lm,lp);
      outVars["llcsip"]=llcsi.first;
      outVars["llcsim"]=llcsi.second;
      outVars["xpt"]=X.Pt();
      outVars["xeta"]=X.Eta();
      outVars["xphi"]=X.Phi(); 
      outVars["mx"]=X.M();
      outVars["xid"]=xid;
      outVars["xht"]=xht;
      outVars["mtx"]=mtx;

      outVars["fpt"]=F.Pt();
      outVars["feta"]=F.Eta();
      outVars["fphi"]=F.Phi();
      outVars["mf"]=F.M();
      outVars["fht"]=fht;      
      outVars["facopl"]=facopl;
      outVars["fcosthetaStarCS"]=fcosthetaStarCS;
      outVars["fphiStar"]=fphiStar;
      outVars["fMR"]=fMR;
      outVars["fR"]=fR;
      std::pair<Float_t,Float_t> fcsi=calcCsi(dil,X);
      outVars["fcsip"]=fcsi.first;
      outVars["fcsim"]=fcsi.second;

      outVars["nb"]=bJets.size();
      outVars["nj"]=lightJets.size();
      outVars["nl"]=leptons.size();
      outVars["ng"]=photons.size();
      outVars["nch"]=nch;
      outVars["ht"]=scalarht;
      outVars["htb"]=scalarhtb;
      outVars["htj"]=scalarhtj;
      outVars["closestkdz"]=closestTrackDZ;
      outVars["puppirecoil_pt"]=puppiRecoil.Mod();
      outVars["puppirecoil_phi"]=puppiRecoil.Phi();
      outVars["puppirecoil_spher"]=fabs(puppiRecoil.Mod())/puppiRecoilHt;

      nRPtk=0;
      if (ev.isData) {
        const edm::EventID ev_id( ev.run, ev.lumi, ev.event );        
        const ctpps::conditions_t lhc_cond = lhc_conds.get( ev_id );
        const double xangle = lhc_cond.crossing_angle;
        for (int ift=0; ift<ev.nfwdtrk; ift++) {
          // Pot ID consist of arm number, station number, pot number
          const unsigned short pot_raw_id = 100*ev.fwdtrk_arm[ift]+10*ev.fwdtrk_station[ift]+ev.fwdtrk_pot[ift];

          // No alignment parameters or dispersion values available for Pot numbers different than 003,103,023 or 123, see https://github.com/forthommel/CTPPSAnalysisTools/blob/83779a55503dcc377fd994719d4aa42b6f9654e3/data/2017/alignments_30jan2017.txt 
          if (pot_raw_id!=3 && pot_raw_id!=23 && pot_raw_id!=103 && pot_raw_id!=123) continue;

          const ctpps::alignment_t align = ctpps_aligns.get( ev_id, pot_raw_id );
          double xi, xi_error;

          // No dispersion values available for xangles different than 120, 130 or 140 murad. See https://github.com/forthommel/CTPPSAnalysisTools/blob/83779a55503dcc377fd994719d4aa42b6f9654e3/data/2017/dispersions.txt
          if (xangle!=120 && xangle!=130 && xangle!=140) continue; 
          proton_reco.reconstruct(xangle, pot_raw_id, ev.fwdtrk_x[ift]*100+align.x_align, xi, xi_error);

          RPid[nRPtk]=pot_raw_id;
          RPcsi[nRPtk]=xi;
          RPcsiUnc[nRPtk]=xi_error;
          RPxang[nRPtk]=xangle;
          nRPtk++;
        }
      }

      outT->Fill();
    }
    
  
  //close input file
  f->Close();
  
  //save histos to file  
  fOut->cd();
  for (auto& it : ht.getPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  outT->Write();
  fOut->Close();
}
