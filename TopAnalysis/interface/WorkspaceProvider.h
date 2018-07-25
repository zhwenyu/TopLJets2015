//////////////////////////////////////
// Current;y the simplest case:     //
//    - Only signal region          //
//    - Hence no transfer factor    //
//////////////////////////////////////

#ifndef WSPROVIDER_H
#define WSPROVIDER_H
#include <iostream>
#include <sstream>
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooBernstein.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <vector>
#include <sstream>
#include <iostream>
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TColorWheel.h"
#include "TColorGradient.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooRandom.h"
#include <algorithm>
#include "RooWorkspace.h"
#include <fstream>
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "RooRandom.h"
#include <iostream>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TArrow.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "TChain.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooArgList.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooVoigtian.h"
#include "RooPlot.h"
#include "RooArgusBG.h"
#include <stdio.h>
#include <math.h>
#include <map>

using namespace std;
using namespace RooFit;

typedef std::map<TString, std::pair<double, std::pair<double, double> > > MeanErr;
typedef std::map<TString, std::pair <double, double> > YieldsErr;

class WorkspaceProvider{
 public:
 WorkspaceProvider(TString histName, TString channel, TString v, int nbin): hist(histName),chan(channel), boson(v),nBin(nbin){
    this->createHistograms();
    ws      = new RooWorkspace("ws"+chan+boson, "ws"+chan+boson);
    var     = new RooRealVar("var","var",-2,3); // to be tuned later
    ws->import(*var);
  }
  void ProvideWS(){
    RooDataHist * dhData   = new RooDataHist("Data"+chan+boson,"Data"+chan+boson,*var,Import(*hData));
    RooDataHist * dhBkg    = new RooDataHist("Background"+chan+boson,"Background"+chan+boson,*var,Import(*hBkg));
    RooDataHist * dhSignal = new RooDataHist("Signal"+chan+boson,"Signal"+chan+boson,*var,Import(*hSig));
    ws->import(*dhData);
    ws->import(*dhBkg);
    ws->import(*dhSignal);
    stringstream s;
    s << hBkg->Integral();
    ws->factory("Background"+chan+boson+"_norm["+s.str()+"]");
    s.str("");
    s << hSig->Integral();
    ws->factory("Signal"+chan+boson+"_norm["+s.str()+"]");   

    TFile * fOut = new TFile("Channel_"+chan+boson+".root","recreate");
    fOut->cd();
    ws->Write();
    fOut->Save();
    fOut->Close();
  }

  void makeCard(YieldsErr YieldErrors, double sigEff=1, double bkgEff=1){
    // No shape uncertainty yet!
    TString outname = chan+".txt";
    ofstream myfile;
    myfile.setf(ios_base::fixed);
    myfile.precision(4);
    myfile.open(outname);

    myfile << "imax 1  number of categories" << endl;
    myfile << "jmax 1  number of samples minus 1" << endl;
    myfile << "kmax *  number of nuisance parameters (sources of systematical uncertainties)" << endl;

    myfile << "\n------------" << endl;
    myfile << "shapes\tSignal\t"<<chan<<boson<<"\tChannel_"<<chan<<boson<<".root ws"<<chan<<boson<<":Signal"<<chan<<boson << endl; 
    myfile << "shapes\tBkg\t"   <<chan<<boson<<"\tChannel_"<<chan<<boson<<".root ws"<<chan<<boson<<":Background"<<chan<<boson << endl;
    myfile << "shapes\tdata_obs\t"<<chan<<boson<<"\tChannel_"<<chan<<boson<<".root ws"<<chan<<boson<<":Data"<<chan<<boson << endl;
    myfile << "------------" << endl;
    myfile << "bin\t"<<chan<<boson << endl;
    myfile << "observation\t-1.0" << endl;
    myfile << "------------" << endl;

    myfile << "bin\t"<<chan<<boson<<"\t"<<chan<<boson << endl;
    myfile << "process\tSignal\tBkg" << endl;
    myfile << "process\t0\t1" << endl;
    myfile << "rate\t"<<hSig->Integral()<<"\t"<<hBkg->Integral()<< endl;
    myfile << "------------" << endl;
    
    for (auto& x : YieldErrors) {
      if(x.second.first == x.second.second )
	myfile << x.first << "\tlnN\t" << x.second.first << "\t-" << endl;
      else
	myfile << x.first << "\tlnN\t" << x.second.first << "/" << x.second.second << "\t-" << endl;
    }
    myfile.close();
  }

  void createHistograms(){
    TString bkgs[]={"Top+VV","DY","QCD","#gamma+jets"};
    TString signal = "EWK #gammaJJ";
    if(boson == "MM")
      signal = "EWK ZJJ";
    TDirectory * dir = (TDirectory*)((TFile::Open("plotter_"+chan+".root"))->Get(chan+boson+"_"+hist));
    dir->ls();
    
    hBkg = (TH1F*)dir->Get(chan+boson+"_"+hist+"_"+bkgs[0]);
    hBkg->SetNameTitle("Background","Background");
    int nRebin =(int)((double)hBkg->GetXaxis()->GetNbins()/(double)nBin);
    for(int i = 1; i < 4; i++){
      TH1F * tmp = (TH1F*)dir->Get(chan+boson+"_"+hist+"_"+bkgs[i]);
      hBkg->Add(tmp);
    }
    hBkg->Rebin(nRebin);
    
    hSig = (TH1F*)dir->Get(chan+boson+"_"+hist+"_"+signal);
    hSig->SetNameTitle("Signal","Signal");
    hSig->Rebin(nRebin);
    
    hData = (TH1F*)dir->Get(chan+boson+"_"+hist);
    hData->SetNameTitle("Data","Data");
    hData->Rebin(nRebin);
    
    TFile* fBkg = new TFile("Background_"+chan+boson+".root","recreate");
    fBkg->cd();
    hBkg->Write();
    fBkg->Close();
    
    TFile* fSig = new TFile("Signal_"+chan+boson+".root","recreate");
    fSig->cd();
    hSig->Write();
    fSig->Close();
    
    TFile* fData = new TFile("Data_"+chan+boson+".root","recreate");
    fData->cd();
    hData->Write();
    fData->Close();
  }

  void creatTFHists(){
    TString bkgs[]={"Top+VV","DY","QCD","#gamma+jets"};
    TString signal = "EWK #gammaJJ";
    if(boson == "MM")
      signal = "EWK ZJJ";
    TDirectory * dirMM = (TDirectory*)((TFile::Open("plotter_"+chan+".root"))->Get(chan+"MM_"+hist));
    TDirectory * dirA  = (TDirectory*)((TFile::Open("plotter_"+chan+".root"))->Get(chan+"A_"+hist));

    bkgTF = (TH1F*)dirA->Get(chan+"A_"+hist+"_#gamma+jets");
    bkgTF->SetNameTitle("BackgroundTF","BackgroundTF");
    TH1F * tmp = (TH1F*)dirMM->Get(chan+"MM_"+hist+"_DY");
    bkgTF->Sumw2();
    tmp->Sumw2();
    int nRebin =(int)((double)bkgTF->GetXaxis()->GetNbins()/(double)nBin);
    bkgTF->Rebin(nRebin);
    tmp->Rebin(nRebin);
    bkgTF->Divide(tmp);
    
    sigTF = (TH1F*)dirA->Get(chan+"A_"+hist+"_EWK #gammaJJ");
    sigTF->SetNameTitle("SignalTF","SignalTF");
    tmp = (TH1F*)dirMM->Get(chan+"MM_"+hist+"_EWK ZJJ");
    sigTF->Sumw2();
    tmp->Sumw2();
    sigTF->Rebin(nRebin);
    tmp->Rebin(nRebin);
    sigTF->Divide(tmp);
    
    TFile* fTF = new TFile("TransferFactors_"+chan+".root","recreate");
    fTF->cd();
    sigTF->Write();
    bkgTF->Write();
    fTF->Close();    
    delete tmp;
  }
 private:
  TString hist,chan,boson;
  TH1F * hBkg, * hSig, * hData;
  TH1F * sigTF, * bkgTF;
  RooRealVar * var;
  RooWorkspace * ws;
  int nBin;
};

YieldsErr splitter(const std::string &s, char delim1 , char delim2){
  YieldsErr ret;
  int nSyst = count(s.begin(),s.end(),delim1) + 1;
  string mystr = s;
  for(int i = 0; i < nSyst; i++){
    int pos1 = mystr.find(delim1);
    string syst = mystr.substr(0,pos1-1);
    string name = syst.substr(0,syst.find(delim2));
    syst = syst.substr(syst.find(delim2)+1);
    double down = (double)atof(syst.substr(0,syst.find(delim2)).c_str());
    syst = syst.substr(syst.find(delim2)+1);
    double up = (double)atof(syst.c_str());
    ret[name]=make_pair(down,up);
  }
  cout << ret.size()<<endl;
  return ret;
}
#endif
