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


//
void RunTTZAnalyzer(TString filename,
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
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
 
  //READ TREE FROM FILE
  MiniEvent_t ev;
  MiniEvent_t& evref = ev;
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

  ht.addHist("genZm", new TH1F("genZmb",";M(l,l');Events",50,60,100));
  ht.addHist("nbGenLeptons", new TH1F("nbGenLeptons",";nb. gen. leptons;Events",10,-.5,9.5));

  std::cout << "init done" << std::endl;

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, false, triggerList);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {
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
           
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      TString chTag = selector.flagFinalState(ev);
      if(chTag=="") continue;
      std::vector<Particle> &leptons     = selector.getSelLeptons(); 
      std::vector<Jet>      &jets        = selector.getJets();

      //select a Z candidate out of the leptons
      std::vector<std::pair<int,int> > zCandidates;
      for(size_t il=0; il<leptons.size(); il++)
        {
          for(size_t jl=il+1; jl<leptons.size(); jl++)
            {
              if( abs(leptons[il].id())!=abs(leptons[jl].id()) ) continue;
              if( leptons[il].charge()*leptons[jl].charge()>=0) continue;
              TLorentzVector dil=leptons[il].p4()+leptons[jl].p4();
              if( fabs(dil.M()-91)>20 ) continue;
              zCandidates.push_back( std::pair<int,int>(il,jl) );
            }
        }

     
      //separate jets according to flavour
      std::vector<Jet> bJets,lightJets;
      for(size_t ij=0; ij<jets.size(); ij++) 
        {
          if(jets[ij].flavor()==5) bJets.push_back(jets[ij]);
          else                     lightJets.push_back(jets[ij]);
        }

      //Dummy weights
      float wgt1(1.0);
      std::vector<double>plotwgts1(1,wgt1);

      ////////////////
      //EVENTS WGTS//
      //////////////
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
        std::vector<Particle> trigLepton(1,leptons[0]);
        EffCorrection_t trigSF = lepEffH.getTriggerCorrection(trigLepton, period);
        EffCorrection_t  selSF= lepEffH.getOfflineCorrection(leptons[0], period);

        wgt *= puWgt*trigSF.first*selSF.first;
      
        // generator level weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);

        //update weight for plotter
        plotwgts[0]=wgt;
      }

      ////////////////////
      //EVENT SELECTION//
      //////////////////
      ht.fill("cutsBefore",.5,plotwgts1);

      if(leptons.size()!=3){

	ht.fill("cuts",.5,plotwgts1);
	ht.fill("cutsPlotwgts",.5,plotwgts);
	continue;}

      ht.fill("cutsBefore",1.5,plotwgts1);

      if(zCandidates.size()!=1) {

	ht.fill("cuts",1.5,plotwgts1);
	ht.fill("cutsPlotwgts",1.5,plotwgts);
	continue;}

      bool passJets(lightJets.size()>=2);
      bool passBJets(bJets.size()>=2);

      ht.fill("cutsBefore",2.5,plotwgts1);

      if(!passJets || !passBJets) {

	ht.fill("cuts",2.5,plotwgts1);
	ht.fill("cutsPlotwgts",2.5,plotwgts);
	continue;}
  
      //Find index of 3rd lepton
      int xlep;
      if((zCandidates[0].first + zCandidates[0].second) == 1) xlep = 2;
      else if((zCandidates[0].first + zCandidates[0].second) == 2) xlep = 1;
      else xlep = 0;
    
      //Select 2 highest pT b jets
      std::vector<std::pair<int,int> > bJetsMax;
      double pTmax = bJets[0].pt()+bJets[1].pt();

      int bindex1 = 0;
      int bindex2 = 1;

      //Find the indexs of the two b jets whose pT sum is maximum
      for(size_t k=1; k<bJets.size(); k++){
	for(size_t j=k+1; j<bJets.size(); j++){

	  if((bJets[k].pt()+bJets[j].pt())>=pTmax){

	    pTmax = bJets[k].pt()+bJets[j].pt();
	    bindex1=k;
	    bindex2=j;

	  }

	}
      }

      bJetsMax.push_back(std::pair<int,int>(bindex1,bindex2));

      //Find the 2 light jets that reconstruct the W mass
      std::vector<std::pair<int,int> > Wcandidate;
      double Wmin = fabs((jets[0].p4()+jets[1].p4()).M()-80.4);

      int index1=0;
      int index2=1;

      for(size_t l=0; l<jets.size(); l++){
	for(size_t n=l+1; n<jets.size(); n++){

	  TLorentzVector invMjets = jets[l].p4()+jets[n].p4();
	  ht.fill("mjjTotal",invMjets.M(),plotwgts);

	  if( fabs(invMjets.M()-80.4)<Wmin ){

	    Wmin = fabs( invMjets.M()-80.4 );
	    index1=l;
	    index2=n;
	    
	  }

	}

      }

      Wcandidate.push_back(std::pair<int,int>(index1,index2));

      //Reconstruct Z kinematics from 2 leptons
      Particle zlepton1=leptons[ zCandidates[0].first ];
      Particle zlepton2=leptons[ zCandidates[0].second ];
      TLorentzVector z=zlepton1.p4()+zlepton2.p4();
 
      //W kinematics from light jets
      Jet Wjet1 = jets[ Wcandidate[0].first ];
      Jet Wjet2 = jets[ Wcandidate[0].second ];
      TLorentzVector w = Wjet1.p4()+Wjet2.p4();

      //Reconstruct neutrino kinematics
      TLorentzVector met(0.,0.,0.,0.);
      met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
      neutrinoPzComputer.SetMET(met);
      neutrinoPzComputer.SetLepton(leptons[xlep].p4());
      float neutrinoPz = neutrinoPzComputer.Calculate();
      TLorentzVector nu(met.Px(),met.Py(),neutrinoPz,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(neutrinoPz,2)));

      //W kinematics from lepton and neutrino
      TLorentzVector xlepton = leptons[xlep].p4();
      TLorentzVector w2 = xlepton+nu;

      //Reconstruct ttbar kinematics from all previous objects
      Jet b1 = bJets[bJetsMax[0].first];
      Jet b2 = bJets[bJetsMax[0].second];
      Jet j1 = jets[index1];
      Jet j2 = jets[index2];
      Particle nuLep = leptons[xlep];
      TLorentzVector ttbar = j1.p4()+j2.p4()+b1.p4()+b2.p4()+nu+nuLep.p4();

      TLorentzVector ttbarZ = z+ttbar;

      ht.fill("ttbarPt",ttbar.Pt(),plotwgts);
      ht.fill("ttbarEta",ttbar.PseudoRapidity(),plotwgts);
      ht.fill("ttbarPhi",ttbar.Phi(),plotwgts);
      ht.fill("ttbarM",ttbar.M(),plotwgts);
      ht.fill("ttbar+Z",ttbarZ.M(),plotwgts);
   
      /////////////////////
      //GEN EVENTS LEVEL//
      ///////////////////
      std::vector<Particle> genLeptonsB = selector.getGenLeptons(ev,20.,2.5);
      std::vector<Particle> &genLeptons = genLeptonsB;
      //Find Z candidate from generated leptons
      std::vector<std::pair<double,double> > zGenCandidates;
      zGenCandidates.push_back(std::pair<double,double>(-1,-1));
      for(size_t il=0; il<genLeptons.size(); il++)
        {
          for(size_t jl=il+1; jl<genLeptons.size(); jl++)
            {
              if( abs(genLeptons[il].id())!=abs(genLeptons[jl].id()) ) continue;
              if( genLeptons[il].charge()*genLeptons[jl].charge()>=0) continue;
              TLorentzVector dgenll=genLeptons[il].p4()+genLeptons[jl].p4();
              if( fabs(dgenll.M()-91)>20 ) continue;
              zGenCandidates[0].first = il;
	      zGenCandidates[0].second = jl;
            }
	    }

      std::cout<<zGenCandidates[0].first<<endl;

      //Reconstruct genZ kinematics
      /* Particle genZlepton1 = genLeptons[ zGenCandidates[0].first  ];
      Particle genZlepton2 = genLeptons[ zGenCandidates[0].second  ];
      TLorentzVector genZ = genZlepton1.p4()+genZlepton2.p4();*/

      //control histograms
      ht.fill("nvtx",     ev.nvtx,        plotwgts);
      ht.fill("nbjets", bJets.size(), plotwgts);
      ht.fill("njets",  jets.size(),  plotwgts);

     
      ht.fill("mll",     z.M(),                                  plotwgts);
      ht.fill("ptll",    z.Pt(),                                 plotwgts);
      ht.fill("dphill",  zlepton1.p4().DeltaPhi(zlepton2.p4()),  plotwgts);

      //Kinematics of each lepton associated with the Z
      ht.fill("zlepton1Pt",zlepton1.pt(),plotwgts);
      ht.fill("zlepton1Eta",zlepton1.eta(),plotwgts);
      ht.fill("zlepton1Phi",zlepton1.phi(),plotwgts);

      ht.fill("zlepton2Pt",zlepton2.pt(),plotwgts);
      ht.fill("zlepton2Eta",zlepton2.eta(),plotwgts);
      ht.fill("zleptons2Phi",zlepton2.phi(),plotwgts);

      //Number of leptons in the event
      ht.fill("nbLeptons",leptons.size(),plotwgts);

      //Kinematics of the light jets
      for(size_t ji=0; ji<jets.size();ji++){

	ht.fill("jetsPt",jets[ji].pt(),plotwgts);
	ht.fill("jetsEta",jets[ji].eta(),plotwgts);
}

      //Kinematics of the b jets
      for(size_t jl=0; jl<bJets.size();jl++){

	ht.fill("bJetsPt",bJets[jl].pt(),plotwgts);
	ht.fill("bJetsEta",bJets[jl].eta(),plotwgts);


}

      ht.fill("mjj", w.M(), plotwgts);

      ht.fill("mlnu", w2.M(), plotwgts);

      //Gen events fill histograms
      ht.fill("nbGenLeptons",genLeptons.size(),plotwgts);
      //ht.fill("genZm",genZ.M(),plotwgts);

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

