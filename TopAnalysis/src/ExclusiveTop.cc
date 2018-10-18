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
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

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
//add jet scale uncertainty
//PPS json
//lumi, puweighting, genweighting
//xsec in the analysis json

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
  
  MiniEvent_t ev;

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

  bool hasETrigger,hasMTrigger,hasMMTrigger,hasEETrigger,hasEMTrigger;
  outT->Branch("hasETrigger",&hasETrigger,"hasETrigger/O");
  outT->Branch("hasMTrigger",&hasMTrigger,"hasMTrigger/O");
  outT->Branch("hasEMTrigger",&hasEMTrigger,"hasEMTrigger/O");
  outT->Branch("hasMMTrigger",&hasMMTrigger,"hasMMTrigger/O");
  outT->Branch("hasEETrigger",&hasEETrigger,"hasEETrigger/O");

  bool isSS,isSF,isZ;
  outT->Branch("isSS",&isSS,"isSS/O");
  outT->Branch("isSF",&isSF,"isSF/O");
  outT->Branch("isZ",&isZ,"isZ/O");
  
  ADDVAR(&ev.rho,"rho","F",outT);
  ADDVAR(&ev.pf_closestDZnonAssoc,"closestDZnonAssoc","F",outT);
  ADDVAR(&ev.pf_ch_wgtSum,"nch","F",outT);
  ADDVAR(&ev.met_pt[1],"met_pt","F",outT);
  ADDVAR(&ev.met_phi[1],"met_phi","F",outT);
  ADDVAR(&ev.met_sig[1],"met_sig","F",outT);
  TString fvars[]={"evwgt", "dilcode", 
                   "l1pt", "l1eta", "l1phi", "ml1", "l1id", "mt1",
                   "l2pt", "l2eta", "l2phi", "ml2", "l2id", "mt2",
                   "llpt", "lleta", "llphi", "mll", "llht", "llacopl", "llcosthetaCS", "llphistar", "llMR", "llR", "mtll", 
                   "llcsip", "llcsipUnc", "llcsim", "llcsimUnc",
                   "xpt",  "xeta",  "xphi",  "mx",  "xid", "xht", "mtx",
                   "fpt",  "feta",  "fphi",  "mf",  "fht", "facopl", "fcosthetaCS", "fphistar", "fMR","fR", 
                   "fcsip", "fcsipUnc", "fcsim", "fcsimUnc",
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
    
  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;

  //LEPTON EFFICIENCIES
  EfficiencyScaleFactorsWrapper lepEffH(filename.Contains("Data13TeV"),era);

  //B-TAG CALIBRATION
  BTagSFUtil btvSF(era,"DeepCSV",BTagEntry::OperatingPoint::OP_MEDIUM,"",0);
  
   //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("puwgtctr",     new TH1F("puwgtctr",    ";Weight sums;Events",2,0,2));
  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",50,-0.5,49.5));
  ht.addHist("nlep",         new TH1F("nlep",        ";Lepton multipliciy;Events",3,2,5));
  ht.addHist("nljets",       new TH1F("nljets",      ";light jet multiplicity;Events",6,0,6)); 
  ht.addHist("nbjets",       new TH1F("nbjets",      ";b jet multiplicity;Events",5,0,5));
  ht.addHist("lmpt",         new TH1F("lmpt",        ";Lepton 1 transverse momentum [GeV]",50,20,200));
  ht.addHist("lmeta",        new TH1F("lmeta",       ";Lepton 1 pseudo-rapidity",10,0,2.5));
  ht.addHist("lppt",         new TH1F("lppt",        ";Lepton 2 transverse momentum [GeV]",50,20,200));
  ht.addHist("lpeta",        new TH1F("lpeta",       ";Lepton 2 pseudo-rapidity",10,0,2.5));
  ht.addHist("mll",          new TH1F("mll",         ";Dilepton invariant mass [GeV]",50,20,200));
  ht.addHist("ptll",         new TH1F("ptll",        ";Dilepton transverse momentum [GeV]",50,0,200));  
  ht.addHist("phistar",      new TH1F("phistar",     ";log_{10}(dilepton #phi^{*})",50,-3,1));
  ht.addHist("costhetaCS",   new TH1F("costhetaCS",  ";Dilepton cos#theta^{*}_{CS}",50,-1,2));
  ht.addHist("acopl",        new TH1F("acopl",       ";Acoplanarity",50,0,1.0));

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
      hasETrigger=(selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v", ev.triggerBits));
      hasMTrigger=(selector.hasTriggerBit("HLT_IsoMu24_v",     ev.triggerBits) ||
                   selector.hasTriggerBit("HLT_IsoMu24_2p1_v", ev.triggerBits) ||
                   selector.hasTriggerBit("HLT_IsoMu27_v",     ev.triggerBits) );
      hasMMTrigger=(selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                  ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",          ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",        ev.triggerBits) );
      hasEETrigger=(selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",             ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          ev.triggerBits) );
      hasEMTrigger=(selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev.triggerBits) ||
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
      std::vector<Particle> leptons        = selector.selLeptons(flaggedLeptons,SelectionTool::TIGHT,SelectionTool::TIGHT);
      if(leptons.size()<2) continue;
      if(leptons[0].Pt()<30 || fabs(leptons[0].Eta())>2.1) continue;

      //Z finder (leading pT leptons)
      int l1idx(0),l2idx(1);
      for(size_t il=0; il<leptons.size(); il++)
        for(size_t jl=1; jl<leptons.size(); jl++) {
          if(abs(leptons[il].id())!=abs(leptons[jl].id())) continue;
          float mll((leptons[il]+leptons[jl]).M());
          if( fabs(mll-91)>10 ) continue;
          l1idx=il;
          l2idx=jl;
          break;
        }
      isSF=( leptons[l1idx].id()==leptons[l2idx].id() );
      isSS=( leptons[l1idx].charge()*leptons[l2idx].charge() > 0 );

      TLorentzVector lm(leptons[l1idx].charge()>0 ? leptons[l1idx] : leptons[l2idx]);
      float lmScaleUnc(leptons[l1idx].charge()>0 ? leptons[l1idx].scaleUnc() : leptons[l2idx].scaleUnc());
      TLorentzVector lp(leptons[l1idx].charge()>0 ? leptons[l2idx] : leptons[l1idx]);
      float lpScaleUnc(leptons[l1idx].charge()>0 ? leptons[l2idx].scaleUnc() : leptons[l1idx].scaleUnc());
      if(isSS)  { 
        lm=leptons[l1idx]; 
        lmScaleUnc=leptons[l1idx].scaleUnc();
        lp=leptons[l2idx];
        lpScaleUnc=leptons[l2idx].scaleUnc();
      }
      TLorentzVector dil(lm+lp);
      Float_t dilScaleUnc=TMath::Sqrt( pow(lm.Pt()*lmScaleUnc,2)+pow(lp.Pt()*lpScaleUnc,2) )/dil.Pt();
      float mll(dil.M());
      if(mll<20) continue;
      isZ=( isSF && !isSS && fabs(mll-91)<10);
 
      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt[1],0,ev.met_phi[1],0.);
      Float_t metScaleUnc(1./ev.met_sig[1]);

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
          float deepcsv(ev.j_deepcsv[idx]);
          int jid=ev.j_id[idx];
          bool passLoosePu((jid>>2)&0x1);
          if(!passLoosePu) continue;
          if(deepcsv>0.4941) { bJets.push_back(allJets[ij]);     scalarhtb+=allJets[ij].pt();  }
          else                             { lightJets.push_back(allJets[ij]); scalarhtj+= allJets[ij].pt(); }
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

      Int_t dilcode=leptons[l1idx].id()*leptons[l2idx].id();

      TString dilCat(isSS ? "ss" : "os" );
      dilCat+=!isSF ? "em" : dilcode==121 ? "ee" : "mm";
      if(isZ) dilCat += "z";
      cats.push_back(dilCat);       
      
      TLorentzVector X(met);
      float xScaleUnc(metScaleUnc);
      float xEtaUnc(0);
      Int_t xid(0);
      Float_t xht(0);
      Float_t mt3l(-1);      
      cats.push_back(dilCat+"inv");
      if(leptons.size()>2){
        cats[2]=dilCat+"3l"; 
        int l3idx(1);
        if(l1idx==0) {
          if(l2idx==1) l3idx=2;
        }
        else if(l1idx==1) {
          l3idx=0;
        }

        neutrinoPzComputer.SetMET(met);
        neutrinoPzComputer.SetLepton(leptons[l3idx].p4());
        float nupz=neutrinoPzComputer.Calculate();
        float metShiftUp(fabs(1+metScaleUnc));
        TLorentzVector upMET(metShiftUp*met);
        neutrinoPzComputer.SetMET(upMET);
        float nupzUp=neutrinoPzComputer.Calculate();
        float metShiftDn(fabs(1-metScaleUnc));
        TLorentzVector dnMET(metShiftDn*met);
        neutrinoPzComputer.SetMET(dnMET);        
        float nupzDn=neutrinoPzComputer.Calculate();

        TLorentzVector neutrinoP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));
        TLorentzVector neutrinoP4Up(upMET.Px(),upMET.Py(),nupzUp ,TMath::Sqrt(TMath::Power(upMET.Pt(),2)+TMath::Power(nupzUp,2)));
        TLorentzVector neutrinoP4Dn(dnMET.Px(),dnMET.Py(),nupzDn ,TMath::Sqrt(TMath::Power(dnMET.Pt(),2)+TMath::Power(nupzDn,2)));

        X=leptons[l3idx]+neutrinoP4;
        xScaleUnc=TMath::Sqrt(
                              pow(leptons[l3idx].scaleUnc()*leptons[l3idx].Pt(),2)+
                              pow(metScaleUnc*neutrinoP4.Pt(),2)
                              )/X.Pt();
        xEtaUnc=max( fabs((leptons[l3idx]+neutrinoP4Up).Eta()-X.Eta()),
                     fabs((leptons[l3idx]+neutrinoP4Dn).Eta()-X.Eta()));
        xht=X.Pt();
        xid=0;
        mt3l=computeMT(leptons[l3idx],met);
      }
      else if(photons.size()>0) {
        cats[2]=dilCat+"a";
        X=photons[0].p4();
        xScaleUnc=photons[0].scaleUnc();
        xid=22;
        xht=X.Pt();
      }
      else if (jets.size()>2) {
        cats[2]=dilCat+"jj";
        X=jets[0]+jets[1];
        xScaleUnc=TMath::Sqrt(pow(jets[0].getScaleUnc()*jets[0].Pt(),2)+
                              pow(jets[1].getScaleUnc()*jets[1].Pt(),2))/X.Pt();
        xid=2121;
        xht=jets[0].Pt()+jets[1].Pt();
        if(bJets.size()>1) {          
          cats[2]=dilCat+"bb";
          X=bJets[0]+bJets[1];
          xScaleUnc=TMath::Sqrt(pow(bJets[0].getScaleUnc()*bJets[0].Pt(),2)+
                                pow(bJets[1].getScaleUnc()*bJets[1].Pt(),2))/X.Pt();
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
      if(mt3l>0) mtx=mt3l; //for 3-leptons use the MT(3rd lep,MET) instead

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
      outVars["evwgt"]=plotwgts[0];
      outVars["dilcode"]=float(dilcode);

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

      ValueCollection_t llcsi=calcCsi(lm,lmScaleUnc,lp,lpScaleUnc);
      outVars["llcsip"]=llcsi[0].first;
      outVars["llcsipUnc"]=llcsi[0].second;
      outVars["llcsim"]=llcsi[1].first;
      outVars["llcsimUnc"]=llcsi[1].second;
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
      outVars["fcosthetaCS"]=fcosthetaStarCS;
      outVars["fphistar"]=fphiStar;
      outVars["fMR"]=fMR;
      outVars["fR"]=fR;
      ValueCollection_t fcsi=calcCsi(dil,dilScaleUnc,X,xScaleUnc,xEtaUnc);
      outVars["fcsip"]=fcsi[0].first;
      outVars["fcsipUnc"]=fcsi[0].second;
      outVars["fcsim"]=fcsi[1].first;
      outVars["fcsimUnc"]=fcsi[1].second;
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
