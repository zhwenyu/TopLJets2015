#include "TopLJets2015/TopAnalysis/interface/JECTools.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TH2F.h"


//
FactorizedJetCorrector *getFactorizedJetEnergyCorrector(TString baseDir, bool isMC)
{
  gSystem->ExpandPathName(baseDir);
  
  //order matters: L1 -> L2 -> L3 (-> Residuals)
  std::vector<std::string> jetCorFiles;
  TString pf("Summer15_25nsV7_");
  pf += (isMC ? "MC" : "DATA");
  std::cout << "Jet energy corrections from " << baseDir+"/"+pf+"_*_AK4PFchs.txt" << std::endl;
  jetCorFiles.push_back((baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt").Data());
  jetCorFiles.push_back((baseDir+"/"+pf+"_L2Relative_AK4PFchs.txt").Data());
  jetCorFiles.push_back((baseDir+"/"+pf+"_L3Absolute_AK4PFchs.txt").Data());
  if(!isMC) jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
  
  //init the parameters for correction
  std::vector<JetCorrectorParameters> corSteps;
  for(size_t i=0; i<jetCorFiles.size(); i++) corSteps.push_back(JetCorrectorParameters(jetCorFiles[i]));
  
  //return the corrector
  return new FactorizedJetCorrector(corSteps);
}


//Sources :  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt)
{
  std::vector<float> res(3,1.0);

  float ptSF(1.0), ptSF_err(0.0);
  if(TMath::Abs(eta)<0.5)       { ptSF=1.109; ptSF_err = 0.008; }
  else if(TMath::Abs(eta)<0.8)  { ptSF=1.138; ptSF_err = 0.013; }
  else if(TMath::Abs(eta)<1.1)  { ptSF=1.114; ptSF_err = 0.013; }
  else if(TMath::Abs(eta)<1.3)  { ptSF=1.123; ptSF_err = 0.024; }
  else if(TMath::Abs(eta)<1.7)  { ptSF=1.084; ptSF_err = 0.011; }  
  else if(TMath::Abs(eta)<1.9)  { ptSF=1.082; ptSF_err = 0.035; }
  else if(TMath::Abs(eta)<2.1)  { ptSF=1.140; ptSF_err = 0.047; }
  else if(TMath::Abs(eta)<2.3)  { ptSF=1.067; ptSF_err = 0.053; }
  else if(TMath::Abs(eta)<2.5)  { ptSF=1.177; ptSF_err = 0.041; }
  else if(TMath::Abs(eta)<2.8)  { ptSF=1.364; ptSF_err = 0.039; }
  else if(TMath::Abs(eta)<3.0)  { ptSF=1.857; ptSF_err = 0.071; }
  else if(TMath::Abs(eta)<3.2)  { ptSF=1.328; ptSF_err = 0.022; }
  else if(TMath::Abs(eta)<5.0)  { ptSF=1.160; ptSF_err = 0.029; }

  res[0] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF)*(pt-genjpt)))/pt;
  res[1] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF-ptSF_err)*(pt-genjpt)))/pt;
  res[2] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF+ptSF_err)*(pt-genjpt)))/pt;
  
  return res;
}


//apply jet energy resolutions (scaling method)
void smearJetEnergies(MiniEvent_t &ev, std::string option) {
  if(ev.isData) return;
  
  for (int k = 0; k < ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    //smear jet energy resolution for MC
    float genJet_pt(0);
    if(ev.j_g[k]>-1) genJet_pt = ev.g_pt[ ev.j_g[k] ];
    if(genJet_pt>0) {
      smearJetEnergy(jp4,genJet_pt,option);
      ev.j_pt[k]   = jp4.Pt();
      ev.j_eta[k]  = jp4.Eta();
      ev.j_phi[k]  = jp4.Phi();
      ev.j_mass[k] = jp4.M();
    }
  }
}

//apply jet energy resolutions (hybrid method)
void smearJetEnergies(MiniEvent_t &ev, JME::JetResolution* jer, std::string option) {
  if(ev.isData) return;

  TRandom3* random = new TRandom3(0); // random seed
  
  for (int k = 0; k < ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    //smear jet energy resolution for MC
    float genJet_pt(0);
    if(ev.j_g[k]>-1) genJet_pt = ev.g_pt[ ev.j_g[k] ];
    //scaling method for matched jets
    if(genJet_pt>0) {
      smearJetEnergy(jp4,genJet_pt,option);
      ev.j_pt[k]   = jp4.Pt();
      ev.j_eta[k]  = jp4.Eta();
      ev.j_phi[k]  = jp4.Phi();
      ev.j_mass[k] = jp4.M();
    }
    //stochastic smearing for unmatched jets
    else {
      double jet_resolution = jer->getResolution({{JME::Binning::JetPt, ev.j_pt[k]}, {JME::Binning::JetEta, ev.j_eta[k]}, {JME::Binning::Rho, ev.rho}});
      smearJetEnergyStochastic(jp4,random,jet_resolution,option);
      ev.j_pt[k]   = jp4.Pt();
      ev.j_eta[k]  = jp4.Eta();
      ev.j_phi[k]  = jp4.Phi();
      ev.j_mass[k] = jp4.M();
    }
  }
  
  delete random;
}

//
void smearJetEnergy(TLorentzVector &jp4, float genJet_pt,std::string option)
{
  int smearIdx(0);
  if(option=="up") smearIdx=1;
  if(option=="down") smearIdx=2;
  float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt)[smearIdx];
  jp4 *= jerSmear;
}

//
void smearJetEnergyStochastic(TLorentzVector &jp4, TRandom* random, double resolution, std::string option)
{
  int smearIdx(0);
  if(option=="up") smearIdx=1;
  if(option=="down") smearIdx=2;
  float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),0.)[smearIdx];
  float jerFactor = 1 + random->Gaus(0, resolution) * sqrt(std::max(pow(jerSmear, 2) - 1., 0.));
  jp4 *= jerFactor;
}

// See https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources#Main_uncertainties_2016_80X
void applyJetCorrectionUncertainty(MiniEvent_t &ev, JetCorrectionUncertainty *jecUnc, TString jecVar, TString direction) {
  for (int k = 0; k < ev.nj; k++) {
    if ((jecVar == "FlavorPureGluon"  and not (ev.j_flav[k] == 21 or ev.j_flav[k] == 0)) or
        (jecVar == "FlavorPureQuark"  and not (abs(ev.j_flav[k]) <= 3 and abs(ev.j_flav[k]) != 0)) or
        (jecVar == "FlavorPureCharm"  and not (abs(ev.j_flav[k]) == 4)) or
        (jecVar == "FlavorPureBottom" and not (abs(ev.j_flav[k]) == 5)))
      continue;
    
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
    applyJetCorrectionUncertainty(jp4,jecUnc,direction);
    ev.j_pt[k]   = jp4.Pt(); 
    ev.j_eta[k]  = jp4.Eta();
    ev.j_phi[k]  = jp4.Phi();
    ev.j_mass[k] = jp4.M();
  }
}

//
void applyJetCorrectionUncertainty(TLorentzVector &jp4,JetCorrectionUncertainty *jecUnc,TString direction)
{
    jecUnc->setJetPt(jp4.Pt());
    jecUnc->setJetEta(jp4.Eta());
    double scale = 1.;
    if (direction == "up")
      scale += jecUnc->getUncertainty(true);
    else if (direction == "down")
      scale -= jecUnc->getUncertainty(false);
    
    jp4 *= scale;
}

//b fragmentation
std::map<TString, TGraph*> getBFragmentationWeights(TString era) {
  std::map<TString, TGraph*> bfragMap;

  TString bfragWgtUrl(era+"/bfragweights.root");
  gSystem->ExpandPathName(bfragWgtUrl);
  TFile *fIn=TFile::Open(bfragWgtUrl);
  bfragMap["upFrag"] = (TGraph *)fIn->Get("upFrag");
  bfragMap["centralFrag"] = (TGraph *)fIn->Get("centralFrag");
  bfragMap["downFrag"] = (TGraph *)fIn->Get("downFrag");
  bfragMap["PetersonFrag"] = (TGraph *)fIn->Get("PetersonFrag");
  return bfragMap;
}

double computeBFragmentationWeight(MiniEvent_t &ev, TGraph* wgtGr) {
  double weight = 1.;
  for (int i = 0; i < ev.ng; i++) {
    if (abs(ev.g_id[i])==5) weight *= wgtGr->Eval(ev.g_xb[i]);
  }
  return weight;
}

std::map<TString, std::map<int, double> > getSemilepBRWeights(TString era) {
  std::map<TString, TGraph*> bfragMap;
  std::map<TString, std::map<int, double> > brMap;

  TString bfragWgtUrl(era+"/bfragweights.root");
  gSystem->ExpandPathName(bfragWgtUrl);
  TFile *fIn=TFile::Open(bfragWgtUrl);
  bfragMap["semilepbrUp"] = (TGraph *)fIn->Get("semilepbrUp");
  bfragMap["semilepbrDown"] = (TGraph *)fIn->Get("semilepbrDown");
  
  for (auto const& gr : bfragMap) {
    for (int i = 0; i < gr.second->GetN(); ++i) {
      double x,y;
      gr.second->GetPoint(i,x,y);
      brMap[gr.first][x] = y;
    }
  }
  
  return brMap;
}

double computeSemilepBRWeight(MiniEvent_t &ev, std::map<int, double> corr, int pid, bool useabs) {
  double weight = 1.;
  for (int i = 0; i < ev.ng; i++) {
    if (!ev.g_isSemiLepBhad[i]) continue;
    if (corr.count(ev.g_bid[i]) == 0) continue;
    if (!useabs and (pid == 0 or pid == ev.g_bid[i])) weight *= corr[ev.g_bid[i]];
    else if (useabs and (pid == 0 or pid == abs(ev.g_bid[i]))) {
      weight *= (corr[ev.g_bid[i]]+corr[-ev.g_bid[i]])/2.;
    }
  }
  return weight;
}

