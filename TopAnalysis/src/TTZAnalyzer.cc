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
#include "TopLJets2015/TopAnalysis/interface/MttbarAnalyzer.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/TOPJetShape.h"


#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"


using namespace std;

struct TTZEvent{

  bool passSelection;
  //Leptons kinematics
  TLorentzVector l0,zl1,zl2;
  //Jets kinematics
  TLorentzVector wj1,wj2,jb1,jb2;
  //W kinematics
  TLorentzVector wjj,wlnu;
  //ttbar kinematics 
  TLorentzVector ttbar;
  //Z kinematics
  TLorentzVector z;

  //Particle lep;
  std::vector<Particle> Leptons;
  std::vector<Jet> Jets;
  std::vector<Jet> LightJets;
  std::vector<Jet> BJets;

};

TTZEvent selectRecoEvents(SelectionTool &selector,MiniEvent_t &ev,MEzCalculator &neutrinoPzComputer){
  
    TTZEvent reco;
    reco.passSelection=false;
    
    //select a Z candidate out of the leptons
    reco.Leptons = selector.getSelLeptons(); 
    std::vector<std::pair<int,int> > zCandidates;
    for(size_t il=0; il<reco.Leptons.size(); il++)
      {
        for(size_t jl=il+1; jl<reco.Leptons.size(); jl++)
          {
            if( abs(reco.Leptons[il].id())!=abs(reco.Leptons[jl].id()) ) continue;
            if( reco.Leptons[il].charge()*reco.Leptons[jl].charge()>=0) continue;
            TLorentzVector dil=reco.Leptons[il].p4()+reco.Leptons[jl].p4();
            if( fabs(dil.M()-91)>20 ) continue;
            zCandidates.push_back( std::pair<int,int>(il,jl) );
          }
      }
     
    //jet selection
    reco.Jets = selector.getJets();
    for(size_t ij=0; ij<reco.Jets.size(); ij++) 
      {
        if(reco.Jets[ij].flavor()==5) reco.BJets.push_back(reco.Jets[ij]);
        else                     reco.LightJets.push_back(reco.Jets[ij]);
      }
    bool passJets(reco.Jets.size()>=4);
    bool passBJets(reco.BJets.size()>=2);

    //preselect the event =3l, =1Z, >=4j, >=2b
    TString chTag = selector.flagFinalState(ev);
    if(chTag=="") return reco;
    if(reco.Leptons.size()!=3 || zCandidates.size()!=1 || !passJets || !passBJets ) return reco;

    reco.passSelection = true;

    //Reconstruct Z kinematics from 2 leptons
    reco.zl1=reco.Leptons[ zCandidates[0].first ].p4();
    reco.zl2=reco.Leptons[ zCandidates[0].second ].p4();
    reco.z=reco.zl1+reco.zl2;

    //Find index of 3rd lepton
    int xlep;
    if((zCandidates[0].first + zCandidates[0].second) == 1) xlep = 2;
    else if((zCandidates[0].first + zCandidates[0].second) == 2) xlep = 1;
    else xlep = 0;
    //Reconstruct 3rd lepton
    reco.l0 = reco.Leptons[xlep].p4();
    
    //Reconstruct neutrino kinematics
    TLorentzVector met(0.,0.,0.,0.);
    met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
    neutrinoPzComputer.SetMET(met);
    neutrinoPzComputer.SetLepton(reco.l0);
    float neutrinoPz = neutrinoPzComputer.Calculate();
    TLorentzVector nu(met.Px(),met.Py(),neutrinoPz,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(neutrinoPz,2)));
    
    //W kinematics from lepton and neutrino
    reco.wlnu = reco.l0+nu;

    //jets are already ordered by pT, can take the first two
    reco.jb1 = reco.BJets[0].p4();
    reco.jb2 = reco.BJets[1].p4();

    //Find the 2 jets that reconstruct closely the W mass
    float Wmin(9999999.);
    int index1(-1),index2(-1);
    for(size_t l=0; l<reco.Jets.size(); l++){
      for(size_t n=l+1; n<reco.Jets.size(); n++){

        //don't consider the jet if it overlaps with the b-jet candidate
        if( reco.Jets[l].p4().DeltaR( reco.BJets[0].p4() )<0.4 ||
            reco.Jets[l].p4().DeltaR( reco.BJets[1].p4() )<0.4 ||
            reco.Jets[n].p4().DeltaR( reco.BJets[0].p4() )<0.4 ||
            reco.Jets[n].p4().DeltaR( reco.BJets[1].p4() )<0.4 ) continue;
          
        TLorentzVector invMjets = reco.Jets[l].p4()+reco.Jets[n].p4();
        if ( fabs(invMjets.M()-80.4)>Wmin ) continue;
        Wmin = fabs( invMjets.M()-80.4 );
        index1=l;
        index2=n;
      }
    } 
      
    //W kinematics from light jets
    Jet Wjet1 = reco.Jets[index1];
    Jet Wjet2 = reco.Jets[index2];
    reco.wj1 = Wjet1.p4();
    reco.wj2 = Wjet2.p4();
    reco.wjj = Wjet1.p4()+Wjet2.p4();

    //Reconstruct ttbar kinematics from all previous objects
    reco.ttbar = reco.wj1+reco.wj2+reco.jb1+reco.jb2+nu+reco.l0;

    //return the event
    return reco;
}

//
TTZEvent selectGenEvents(SelectionTool &selector,MiniEvent_t &ev){

  TTZEvent gen;
  gen.passSelection=false;

  //Find Z candidate from generated leptons
  gen.Leptons = selector.getGenLeptons(ev,20.,2.5);
  std::vector<std::pair<int,int> > zGenCandidates;
  for(size_t il=0; il<gen.Leptons.size(); il++)
    {
      for(size_t jl=il+1; jl<gen.Leptons.size(); jl++)
        {
          if( abs(gen.Leptons[il].id())!=abs(gen.Leptons[jl].id()) ) continue;
          if( gen.Leptons[il].charge()*gen.Leptons[jl].charge()>=0) continue;
          TLorentzVector dgenll=gen.Leptons[il].p4()+gen.Leptons[jl].p4();
          if( fabs(dgenll.M()-91)>20 ) continue;
          zGenCandidates.push_back(std::pair<int,int>(il,jl));
        }
    }

  //select jets
  gen.Jets    = selector.getGenJets(ev,20.,2.5);
  for(size_t ij=0; ij<gen.Jets.size(); ij++) 
    {
      if(gen.Jets[ij].flavor()==5) gen.BJets.push_back(gen.Jets[ij]);
      else                         gen.LightJets.push_back(gen.Jets[ij]);
    }
  bool passGenJets(gen.Jets.size()>=4);
  bool passGenBJets(gen.BJets.size()>=2);

  //preselect the event
  if(gen.Leptons.size()!=3 || zGenCandidates.size()!=1 || !passGenJets || !passGenBJets ) return gen;

  gen.passSelection = true;

  //Reconstruct Z kinematics from 2 leptons
  gen.zl1=gen.Leptons[ zGenCandidates[0].first ].p4();
  gen.zl2=gen.Leptons[ zGenCandidates[0].second ].p4();
  gen.z=gen.zl1+gen.zl2;
  
  //Find index of 3rd lepton
  int xlep;
  if((zGenCandidates[0].first + zGenCandidates[0].second) == 1) xlep = 2;
  else if((zGenCandidates[0].first + zGenCandidates[0].second) == 2) xlep = 1;
  else xlep = 0;
  gen.l0 = gen.Leptons[xlep].p4();

  //set the b-jets
  gen.jb1=gen.BJets[0].p4();
  gen.jb2=gen.BJets[1].p4();


  //Get top kinematics
  TLorentzVector p4(0.,0.,0.,0.);  
  gen.ttbar = p4;
  for(Int_t igen=0; igen<ev.ngtop; igen++)
    {
      p4.SetPtEtaPhiM(ev.gtop_pt[igen],ev.gtop_eta[igen],ev.gtop_phi[igen],ev.gtop_m[igen]);
      gen.ttbar += p4;
    }
  
  return gen;
}

void RunTTZAnalyzer(TString filename,
                    TString outname,
                    Int_t channelSelection, 
                    Int_t chargeSelection, 
                    TH1F *normH, 
                    TString era,
                    Bool_t debug){
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  TRandom* random = new TRandom(0); // random seed for period selection
  std::vector<RunPeriod_t> runPeriods=getRunPeriods(era);
  bool isTTbar( filename.Contains("_TTJets") or (normH and TString(normH->GetTitle()).Contains("_TTJets")));
  bool isData( filename.Contains("Data") );
  
  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
 
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
  
  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;
 
  //PILEUP WEIGHTING
  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if( !isData ) puWgtGr=getPileupWeightsMap(era,genPU);
    
  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);

  //B-TAG CALIBRATION
  std::map<TString, std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> > btvsfReaders = getBTVcalibrationReadersMap(era, BTagEntry::OP_MEDIUM);
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *>    expBtagEffPy8 = readExpectedBtagEff(era);
  
   //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);

  //Leptons
  ht.addHist("nbLeptons", new TH1F("nbLeptons",";Nb. of leptons;Events",10,-.5,9.5));

  //Leptons associated with the Z
  ht.addHist("zlepton1Pt", new TH1F("zlepton1Pt",";p_{T} [GeV];Events",50,0,500));
  ht.addHist("zlepton1Eta",new TH1F("zlepton1Eta",";eta;Events",50,0,2.5));
  ht.addHist("zlepton1Phi",new TH1F("zlepton1Phi",";phi;Events",50,-3.2,3.2));

  ht.addHist("zlepton2Pt", new TH1F("zlepton2Pt",";p_{T} [GeV];Events",50,0,500));
  ht.addHist("zlepton2Eta",new TH1F("zlepton2Eta",";eta;Events",50,0,2.5));
  ht.addHist("zleptons2Phi",new TH1F("zlepton2Phi",";phi;Events",50,-3.2,3.2));

  //Events lost in each cut
  ht.addHist("cuts", new TH1F("cuts",";cut;Events lost",3,0,3));
  ht.addHist("cutsPlotwgts", new TH1F("cutsPlotwgts",";cuts;Events",3,0,3));
  //Events before each cut
  ht.addHist("cutsBefore", new TH1F("cutsBefore",";cut;Events before cut",3,0,3));
 
  //Jets
  ht.addHist("jetsPt",new TH1F("jetsPt",";p_{T} [GeV];Nb. of jets",50,0,300));
  ht.addHist("jetsEta", new TH1F("jetsEta",";eta;Events",50,0,2.5));
  ht.addHist("bJetsPt",new TH1F("bJetsPt",";p_{T} [GeV];Nb. of b jets",50,0,300));
  ht.addHist("bJetsEta", new TH1F("bJetsEta",";eta;Events",50,0,2.5));

  ht.addHist("puwgtctr",     new TH1F("puwgtctr",    ";Weight sums;Events",2,0,2));
  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",55,-0.5,49.5));
  ht.addHist("njets",        new TH1F("njets",       ";Jet multiplicity;Events",15,-0.5,14.5));
  ht.addHist("nbjets",       new TH1F("nbjets",      ";b jet multiplicity;Events",10,-0.5,9.5));
  ht.addHist("mll",          new TH1F("mll",         ";M(l,l') [GeV];Events",50,91-20,91+20));
  ht.addHist("ptll",         new TH1F("ptll",        ";p_{T}(l,l') [GeV];Events",50,0,250));
  ht.addHist("dphill",       new TH1F("dphill",      ";#Delta#phi(l,l');Events",50,0,3.15));

  //W from light jets
  ht.addHist("mjj", new TH1F("mjj", ";M(j,j') [GeV];Events",50,80.4-20,80.4+20));
  ht.addHist("mjjTotal", new TH1F("mjjTotal", ";M(j,j') [GeV];Events",50,0.,300.));

  //W from lepton+neutrino
  ht.addHist("mlnu", new TH1F("mlnu", ";M(l,nu) [GeV];Events",50,78,85));

  //ttbar
  ht.addHist("ttbarPt", new TH1F("ttbarPt",";p_{T} [GeV];Events",50,0,250));
  ht.addHist("ttbarEta", new TH1F("ttbarEta",";eta;Events",50,0,2.5));
  ht.addHist("ttbarPhi", new TH1F("ttbarPhi",";phi;Events",50,0,3.15));
  ht.addHist("ttbarM", new TH1F("ttbarM",";M(b,b',j,j',l,MET) [GeV]; Events",50,0.,1000.));

  ht.addHist("ttbar+Z", new TH1F("ttbar+Z",";M(ttbar,Z) [GeV];Events",50,300.,1500.));

  ht.addHist("genZm", new TH1F("genZm",";M(l,l');Events",50,60,100));
  ht.addHist("nbGenLeptons", new TH1F("nbGenLeptons",";nb. gen. leptons;Events",10,-.5,9.5));

  ht.addHist("mig_zpt", new TH2F("mig_zpt",";gen;reco;Events",50,0,250,50,0,250));


  std::cout << "init done" << std::endl;

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////

  int l=0;
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, false, triggerList);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {

      l+=1;
      t->GetEntry(iev);
      if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
      
      //assign randomly a run period
      TString period = assignRunPeriod(runPeriods,random);
      
      //////////////////
      // CORRECTIONS //
      ////////////////
      double csvm = 0.8484;
      addBTagDecisions(ev, csvm, csvm);
      if(!ev.isData) smearJetEnergies(ev);

      ////////////////
      //RECO EVENTS//
      //////////////
      TTZEvent reco = selectRecoEvents(selector,ev,neutrinoPzComputer);
      std::vector<Particle> &allLeptons = reco.Leptons;	 

      ///////////////
      //GEN EVENTS//
      /////////////
      TTZEvent gen = selectGenEvents(selector,ev);
      
      //Reconstructed events
      if(reco.passSelection){
      
      //////////////
      //PLOT WGTS//
      ////////////
      float wgt(1.0);
      std::vector<double>plotwgts(1,wgt);
      ht.fill("puwgtctr",0,plotwgts);

	if (!ev.isData) {

	  // norm weight
	  wgt  = (normH? normH->GetBinContent(1) : 1.0);
        
	  // pu weight
	  double puWgt(puWgtGr[period][0]->Eval(ev.g_pu));
	  std::vector<double>puPlotWgts(1,puWgt);
	  ht.fill("puwgtctr",1,puPlotWgts);

	  // lepton trigger*selection weights
	  std::vector<Particle> trigLepton(1,allLeptons[0]);
	  EffCorrection_t trigSF = lepEffH.getTriggerCorrection(trigLepton, period);
	  EffCorrection_t  selSF= lepEffH.getOfflineCorrection(allLeptons[0], period);
	  wgt *= puWgt*trigSF.first*selSF.first;
      
	  // generator level weights
	  wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);
      
	  //update weight for plotter
	  plotwgts[0]=wgt;
	}

	//ttbar kinematics
	ht.fill("ttbarPt",(reco.ttbar).Pt(),plotwgts);
	ht.fill("ttbarEta",(reco.ttbar).PseudoRapidity(),plotwgts);
	ht.fill("ttbarPhi",(reco.ttbar).Phi(),plotwgts);
	ht.fill("ttbarM",(reco.ttbar).M(),plotwgts);
	//ht.fill("ttbar+Z",ttbarZ.M(),plotwgts);
        
	//control histograms
	ht.fill("nvtx",     ev.nvtx,        plotwgts);
	ht.fill("nbjets", (reco.BJets).size(), plotwgts);
	ht.fill("njets",  (reco.Jets).size(),  plotwgts);
	ht.fill("nbLeptons",(reco.Leptons).size(),plotwgts);

	//Z kinematics
	ht.fill("mll",     (reco.z).M(),                                  plotwgts);
	ht.fill("ptll",    (reco.z).Pt(),                                 plotwgts);
	//ht.fill("dphill",  (reco.zl1).Pt().DeltaPhi(reco.zl2),  plotwgts);

	//Kinematics of each lepton associated with the Z
	ht.fill("zlepton1Pt",(reco.zl1).Pt(),plotwgts);
	ht.fill("zlepton1Eta",(reco.zl1).PseudoRapidity(),plotwgts);
	ht.fill("zlepton1Phi",(reco.zl1).Phi(),plotwgts);

	ht.fill("zlepton2Pt",(reco.zl2).Pt(),plotwgts);
	ht.fill("zlepton2Eta",(reco.zl2).PseudoRapidity(),plotwgts);
	ht.fill("zleptons2Phi",(reco.zl2).Phi(),plotwgts);

	/*    //Kinematics of the light jets
	for(size_t ji=0; ji<(reco.Jets).size();ji++){

	  ht.fill("jetsPt",(reco.Jets)[ji].Pt(),plotwgts);
	  ht.fill("jetsEta",(reco.Jets)[ji].PseudoRapidity(),plotwgts);
	}

      //Kinematics of the b jets
	for(size_t jl=0; jl<(reco.BJets).size();jl++){

	  ht.fill("bJetsPt",reco.BJets[jl].Pt(),plotwgts);
	  ht.fill("bJetsEta",reco.BJets[jl].PseudoRapidity(),plotwgts);
	 
	  }*/

      //W kinematics
	ht.fill("mjj", (reco.wjj).M(), plotwgts);
	ht.fill("mlnu", (reco.wlnu).M(), plotwgts);

      }

      //Generator events
      /*   if( reco.passSelection || gen.passSelection ){

	float rec_zpt(reco.z.Pt()),gen_zpt(gen.z.Pt());
	if(reco.passSelection!= true ) rec_zpt=-1;
	if(gen.passSelection!= true )  gen_zpt=-1;
	ht.get2dPlots()["mig_zpt"]->Fill(gen_zpt,rec_zpt,plotwgts);

	}*/

    }

  //close input file
  f->Close();
  //save histos to file  
  fOut->cd();
  for (auto& it : ht.getPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}

