// Author: Stefan Schmitt
// DESY, 14.10.2008

//  Version 17.6, in parallel to changes in TUnfold
//
//  History:
//    Version 17.5, in parallel to changes in TUnfold
//    Version 17.4, in parallel to changes in TUnfold
//    Version 17.3, in parallel to changes in TUnfold
//    Version 17.2, in parallel to changes in TUnfold
//    Version 17.1, in parallel to changes in TUnfold
//    Version 17.0, updated for using the classes TUnfoldDensity, TUnfoldBinning
//    Version 16.1, parallel to changes in TUnfold
//    Version 16.0, parallel to changes in TUnfold
//    Version 15, with automated L-curve scan
//    Version 14, with changes in TUnfoldSys.cxx
//    Version 13, include test of systematic errors
//    Version 12, catch error when defining the input
//    Version 11,  print chi**2 and number of degrees of freedom
//    Version 10,  with bug-fix in TUnfold.cxx
//    Version 9,  with bug-fix in TUnfold.cxx and TUnfold.h
//    Version 8,  with bug-fix in TUnfold.cxx and TUnfold.h
//    Version 7,  with bug-fix in TUnfold.cxx and TUnfold.h
//    Version 6a, fix problem with dynamic array allocation under windows
//    Version 6, bug-fixes in TUnfold.C
//    Version 5, replace main() by testUnfold1()
//    Version 4, with bug-fix in TUnfold.C
//    Version 3, with bug-fix in TUnfold.C
//    Version 2, with changed ScanLcurve() arguments
//    Version 1, remove L curve analysis, use ScanLcurve() method instead
//    Version 0, L curve analysis included here

/*
  This file is part of TUnfold.

  TUnfold is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  TUnfold is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with TUnfold.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <TError.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TF1.h>
#include <TStyle.h>
#include <TVector.h>
#include <TGraph.h>
#include <TFile.h>
#include <TLegend.h>

#include "TUnfoldDensity.h"

// #define VERBOSE_LCURVE_SCAN

using namespace std;

///////////////////////////////////////////////////////////////////////
// 
//  Test program for the classes TUnfold and related
//
//  (1) Generate Monte Carlo and Data events
//      The events consist of
//        signal
//        background
//
//      The signal is a resonance. It is generated with a Breit-Wigner,
//      smeared by a Gaussian
//
//  (2) Unfold the data. The result is:
//      The background level
//      The shape of the resonance, corrected for detector effects
//
//      Systematic errors from the MC shape variation are included
//      and propagated to the result
//
//  (3) fit the unfolded distribution, including the correlation matrix
//
//  (4) save six plots to a file testUnfold1.ps
//        1  2  3
//        4  5  6
//      1: 2d-plot of the matrix decsribing the migrations
//      2: generator-level distributions
//             blue: unfolded data, total errors
//             green: unfolded data, statistical errors
//             red: generated data
//             black: fit to green data points
//      3: detector level distributions
//             blue: unfoldede data, folded back through the matrix
//             black: Monte Carlo (with wrong peal position)
//             blue: data
//      4: global correlation coefficients
//      5: chi**2 as a function of log(tau)
//           the star indicates the final choice of tau
//      6: the L curve
//
///////////////////////////////////////////////////////////////////////

TRandom *rnd=0;

TH2 *gHistInvEMatrix;

TVirtualFitter *gFitter=0;

void divideByBinWidth(TH1* hist) {
  for(Int_t i=1; i<=hist->GetNbinsX(); i++) {
     hist->SetBinContent(i,hist->GetBinContent(i)/hist->GetBinWidth(i));
     hist->SetBinError(i,hist->GetBinError(i)/hist->GetBinWidth(i));
  }
}
void normalize(TH1* hist) {
  hist->Scale(1./hist->Integral());
}

int testUnfold0()
{
  // switch on histogram errors
  TH1::SetDefaultSumw2();

  // show fit result
  gStyle->SetOptFit(1111);
  
  TString basepath = "/afs/cern.ch/work/m/mseidel/TopAnalysis/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/unfolding/fill/";
  
  //============================================
  // generate MC distribution
  //
  TFile file(basepath+"MC13TeV_TTJets.root", "READ");
  TFile datafile(basepath+"Data13TeV.root", "READ");
  
  TString observable  = "mult";
          observable += "_charged";
  TString flavor      = "all";
  bool reg = false;
  
  TH1D *histMgenMC = (TH1D*) file.Get(observable+"_"+flavor+"_responsematrix_px");
  //histMgenMC->Scale(1./histMgenMC->Integral());
  TH1D *histMdetMC = (TH1D*) file.Get(observable+"_"+flavor+"_responsematrix_py");
  //histMdetMC->Scale(1./histMdetMC->Integral());
  TH2D *histMdetGenMC = (TH2D*) file.Get(observable+"_"+flavor+"_responsematrix");
  TH1D* histMdetMCsig = histMdetGenMC->ProjectionY("histMdetNonGenMC", 1, -1);
  TH2D* histMdetGenMCSig = (TH2D*) histMdetGenMC->Clone("histMdetGenMCSig");
  
  double dataMCSF = 36500.*832.;
  histMgenMC->Scale(dataMCSF);
  histMdetMC->Scale(dataMCSF);
  histMdetMCsig->Scale(dataMCSF);
  
  double integral = histMdetGenMC->Integral(0,-1,0,-1);
  std::cout << "integral: " << integral << std::endl;
  double nExpectedEvents = integral*36500*832.;
  std::vector<double> wheel = { 0. };
  for (int i = 0; i < histMdetGenMC->GetNcells()+1; ++i) {
    //std::cout << wheel.back() << std::endl;
    wheel.push_back(wheel.back() + histMdetGenMC->GetBinContent(i));
  }
  
  //histMdetGenMC->Scale(1./histMdetGenMC->Integral());
  //histMdetGenMC->Scale(1000000.);
  
  //Set non-reco bins to 0
  //for (int i = 0; i < histMdetGenMCSig->GetNbinsX()+1; ++i) histMdetGenMCSig->SetBinContent(i, 0, 0.);
  //Set non-gen bins to 0
  for (int i = 0; i < histMdetGenMCSig->GetNbinsY()+2; ++i) histMdetGenMCSig->SetBinContent(0, i, 0.);
  //Set non-reco non-gen bin to 0
  //histMdetGenMC->SetBinContent(0, 0, 36500000.);
  
  TRandom3 random(1);
  
  int nToys = 1;
  
  TH2F* hPull = new TH2F("pull", "pull", histMgenMC->GetNbinsX()+1, 0, histMgenMC->GetNbinsX()+1, 50, -5, 5);
  
//  for (int t = 0; t < nToys; ++t) {
//    std::cout<<"toy number "<<t<<std::endl;
    //============================================
    // generate toy distribution
    //
    //TH2D* histMdetGenToy = (TH2D*) histMdetGenMC->Clone("histMdetGenToy");
    //histMdetGenToy->Scale(dataMCSF);
    //histMdetGenToy->Reset();
    //
    //int drawn = 0;
    //while (drawn < nExpectedEvents) {
    //  double pick = random.Uniform(integral);
    //  //std::cout << pick << ":";
    //  auto low = std::lower_bound(wheel.begin(), wheel.end(), pick);
    //  int position = low - wheel.begin() - 1;
    //  //std::cout << position << " ";
    //  //int position = -1;
    //  //for (unsigned int i = 0; i < wheel.size(); ++i) {
    //  //  if (pick < wheel[i]) break;
    //  //  ++position;
    //  //}
    //  //std::cout << position << " ";
    //  histMdetGenToy->SetBinContent(position, histMdetGenToy->GetBinContent(position) + 1);
    //  ++drawn;
    //}
    //for (int i = 0; i < histMdetGenToy->GetNcells()+1; ++i) {
    //  histMdetGenToy->SetBinError(i, sqrt(histMdetGenToy->GetBinContent(i)));
    //}
    //histMdetGenToy->FillRandom(histMdetGenMC, nExpectedEvents);    
    
    //TH1D *histMdetData = histMdetGenToy->ProjectionY("histMdetData");
    //TH1D *histMgenToy = histMdetGenToy->ProjectionX("histMgenToy");
    
    TH1D *histMdetData = (TH1D*) datafile.Get(observable+"_"+flavor+"_responsematrix_py");

    //double dataMCSF = 36500.; //histMdetData->Integral()/histMdetMC->Integral();
    //dataMCSF = histMdetData->Integral()/histMdetMC->Integral()*100;
    //dataMCSF = 100000.;
    
    //if (not DATA) {
    //  histMdetData->Scale(dataMCSF*832.); // For MC test
    //}
    //histMdetData->Scale(100.);
    //histMdetData->SetBinContent(0, 0.);
    
    for (int i = 0; i < histMdetData->GetNbinsX()+2; ++i) cout << histMdetData->GetBinContent(i) << " +/- " << histMdetData->GetBinError(i) << endl;
    
    //for (int i = 0; i < histMdetNonGenMC->GetNbinsX()+1; ++i) cout << histMdetNonGenMC->GetBinContent(i)*dataMCSF << endl;
    
    // backgrounds
    TH1D *histMdetDataBGSubtracted = (TH1D*) histMdetData->Clone();

    TFile bkgfile(basepath+"MC13TeV_TTJets.root", "READ");
    TH2D *histMdetGenMCbkg = (TH2D*) bkgfile.Get(observable+"_"+flavor+"_responsematrix");
    TH1D *histMdetNonGenMCbkg = histMdetGenMCbkg->ProjectionY("histMdetNonGenMCbkg_MC13TeV_TTJets.root", 0, 0);
    TH1D *histMdetNonGenMCall = histMdetGenMCbkg->ProjectionY("histMdetNonGenMCall_MC13TeV_TTJets.root");
    //double SF = histMdetDataBGSubtracted->Integral(1, -1) * histMdetNonGenMCall->Integral(1,-1);
    //std::cout << "SF=" << SF << std::endl;
    //histMdetDataBGSubtracted->Add(histMdetNonGenMCbkg, -SF*832.);
    //*
    for (int i = 0; i < histMdetDataBGSubtracted->GetNbinsX()+2; ++i) {
      double sf = histMdetNonGenMCbkg->GetBinContent(i) / histMdetNonGenMCall->GetBinContent(i);
      histMdetDataBGSubtracted->SetBinContent(i, (1.-sf)*histMdetDataBGSubtracted->GetBinContent(i));
      histMdetDataBGSubtracted->SetBinError(i, (1.-sf)*histMdetDataBGSubtracted->GetBinError(i));
    }
    //*/
    //*
    std::map<TString, double> backgrounds;
    backgrounds["MC13TeV_SingleTbar_tW"] = 35.85;
    backgrounds["MC13TeV_SingleT_tW"] = 35.85;
    backgrounds["MC13TeV_SingleTbar_t"] = 80.95;
    backgrounds["MC13TeV_SingleT_t"] = 136.02;
    for (auto background : backgrounds) {
      TFile bkgfile(basepath+background.first+".root", "READ");
      TH2D *histMdetGenMCbkg = (TH2D*) bkgfile.Get(observable+"_"+flavor+"_responsematrix");
      TH1D *histMdetNonGenMCbkg = histMdetGenMCbkg->ProjectionY("histMdetNonGenMCbkg_"+background.first, 0, 0);
      histMdetDataBGSubtracted->Add(histMdetNonGenMCbkg, -36500.*background.second);
    }
    //*/
    for (int i = 0; i < histMdetDataBGSubtracted->GetNbinsX()+2; ++i) {
      std::cout << histMdetDataBGSubtracted->GetBinContent(i) << " +/- " << histMdetDataBGSubtracted->GetBinError(i) << std::endl;
    }
    

    
    //histMdetData->SetBinContent(0, 1000000.);


    //=========================================================================
    // set up the unfolding
    // define migration matrix
    TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
    TUnfold::EConstraint constraint = TUnfold::kEConstraintArea;
    TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;
    
    TUnfoldDensity unfold(histMdetGenMCSig,TUnfold::kHistMapOutputHoriz,regMode,constraint,densityFlags);

    // define input and bias scame
    // do not use the bias, because MC peak may be at the wrong place
    // watch out for error codes returned by the SetInput method
    // errors larger or equal 10000 are fatal:
    // the data points specified as input are not sufficient to constrain the
    // unfolding process
    if(unfold.SetInput(histMdetDataBGSubtracted)>=10000) {
      std::cout<<"Unfolding result may be wrong\n";
    }

    //========================================================================
    // the unfolding is done here
    //
    // scan L curve and find best point
    Int_t nScan=100;
    // use automatic L-curve scan: start with taumin=taumax=0.0
    Double_t tauMin=1e-10;
    Double_t tauMax=1e-2;
    TSpline *rhoScan;
    TUnfoldDensity::EScanTauMode tauflag = TUnfoldDensity::kEScanTauRhoAvg;

    int iBest = unfold.ScanTau(nScan,tauMin,tauMax,&rhoScan,tauflag);
    
    // create graphs with one point to visualize best choice of tau
    double tau[1], rho[1];
    
    rhoScan->GetKnot(iBest, tau[0], rho[0]);
    TGraph *bestRho = new TGraph(1,tau,rho);
    double opt_tau = unfold.GetTau();
    cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  TAU,iBest --> " << tau[0] << "  " << opt_tau << "  " << rho[0] << endl;

    //==========================================================================
    // print some results
    //
    //std::cout<<"toy number "<<t<<std::endl;
    std::cout<<"tau="<<unfold.GetTau()<<"\n";
    std::cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
             <<" / "<<unfold.GetNdf()<<"\n";
    std::cout<<"chi**2(sys)="<<unfold.GetChi2Sys()<<"\n";
    
    // Add systematic uncertainties
    double nominalNormalization = histMdetGenMCSig->Integral(1,-1,0,-1) / histMdetGenMCSig->Integral(1,-1,1,-1);
    
    std::vector<TString> uncertainties;
    uncertainties.push_back("fsrup");
    uncertainties.push_back("fsrdn");
    uncertainties.push_back("hdampup");
    uncertainties.push_back("hdampdn");
    uncertainties.push_back("isrup");
    uncertainties.push_back("isrdn");
    uncertainties.push_back("ueup");
    uncertainties.push_back("uedn");
    uncertainties.push_back("herwig");
    
    for (TString uncertainty : uncertainties) {
      TFile uncfile(basepath+"MC13TeV_TTJets_"+uncertainty+".root", "READ");
      TH2D *histMdetGenMCunc = (TH2D*) uncfile.Get(observable+"_"+flavor+"_responsematrix");
      for (int i = 0; i<histMdetGenMCunc->GetNbinsY()+2; ++i) histMdetGenMCunc->SetBinContent(0, i, 0.);
      double normalization = histMdetGenMCunc->Integral(1,-1,0,-1) / histMdetGenMCunc->Integral(1,-1,1,-1);
      //std::cout << "normalization = " << normalization << std::endl;
      for (int i = 0; i<histMdetGenMCunc->GetNbinsX()+2; ++i) histMdetGenMCunc->SetBinContent(i, 0, histMdetGenMCunc->GetBinContent(i)*nominalNormalization/normalization);
      unfold.AddSysError(histMdetGenMCunc, uncertainty, TUnfold::kHistMapOutputHoriz, TUnfoldSys::kSysErrModeMatrix);
    }
    
    // uncertainties from weights
    for (int w = 1; w <= 20; ++w) {
      TH2D *histMdetGenMCunc = (TH2D*) file.Get(observable+"_"+flavor+"_wgt"+std::to_string(w)+"_responsematrix");
      for (int i = 0; i<histMdetGenMCunc->GetNbinsY()+2; ++i) histMdetGenMCunc->SetBinContent(0, i, 0.);
      double normalization = histMdetGenMCunc->Integral(1,-1,0,-1) / histMdetGenMCunc->Integral(1,-1,1,-1);
      //std::cout << "normalization = " << normalization << std::endl;
      for (int i = 0; i<histMdetGenMCunc->GetNbinsX()+2; ++i) histMdetGenMCunc->SetBinContent(i, 0, histMdetGenMCunc->GetBinContent(i)*nominalNormalization/normalization);
      unfold.AddSysError(histMdetGenMCunc, TString(std::to_string(w)), TUnfold::kHistMapOutputHoriz, TUnfoldSys::kSysErrModeMatrix);
    }


    //==========================================================================
    // retreive results into histograms
    
    if (not reg) unfold.DoUnfold(0.);
    else unfold.DoUnfold(opt_tau);
    
    TH1 *histMunfold=unfold.GetOutput("Unfolded");
    
    TH1* histMfoldback = unfold.GetFoldedOutput("FoldedBack");
    
    // get error matrix (input distribution [stat] errors only)
    TH2 *histEmatData=unfold.GetEmatrixInput("EmatData");
    
    TH2 *histEmatStat = unfold.GetEmatrixSysUncorr("EmatStat");

    // get total error matrix:
    //   migration matrix uncorrelated and correlated systematic errors
    //   added in quadrature to the data statistical errors
    TH2 *histEmatTotal=unfold.GetEmatrixTotal("EmatTotal");

    // create data histogram with the total errors
    TH1D *histTotalError=new TH1D("TotalError",";",histMunfold->GetNbinsX(),histMunfold->GetXaxis()->GetXbins()->GetArray());
    histTotalError->SetXTitle(histMunfold->GetXaxis()->GetTitle());
    for(Int_t bin=1;bin<=histMunfold->GetNbinsX();bin++) {
      histTotalError->SetBinContent(bin,histMunfold->GetBinContent(bin));
      histTotalError->SetBinError(bin,TMath::Sqrt(histEmatTotal->GetBinContent(bin,bin)));
      printf("Entries: %9.0g, total error: %6.0f (%5.2f%%), data stat error: %6.0f (%5.2f%%)\n", histTotalError->GetBinContent(bin), histTotalError->GetBinError(bin), histTotalError->GetBinError(bin)/histTotalError->GetBinContent(bin)*100, TMath::Sqrt(histEmatData->GetBinContent(bin,bin)), TMath::Sqrt(histEmatData->GetBinContent(bin,bin))/histTotalError->GetBinContent(bin)*100);
    }


  //=====================================================================
  // plot some histograms
  TCanvas output;
  output.Divide(3,2);

  // Show the matrix which connects input and output
  // There are overflow bins at the bottom, not shown in the plot
  // These contain the background shape.
  // The overflow bins to the left and right contain
  // events which are not reconstructed. These are necessary for proper MC
  // normalisation
  output.cd(1);
  gStyle->SetPaintTextFormat("4.1g");
  histMdetGenMC->Draw("colz,text");

  // draw generator-level distribution:
  //   data (red) [for real data this is not available]
  //   MC input (black) [with completely wrong peak position and shape]
  //   unfolded data (blue)
  output.cd(2);
  
  //normalize(histMgenMC);
  //normalize(histMdetMC);
  //normalize(histMdetData);
  //normalize(histMdetDataBGSubtracted);
  //normalize(histTotalError);
  //normalize(histMunfold);
  
  divideByBinWidth(histMgenMC);
  //histMgenMC->Scale(histMgenMC->Integral()/histMunfold->Integral());
  //histMgenMC->Scale(histMunfold->Integral()/histMgenMC->Integral());
  histMgenMC->GetYaxis()->SetRangeUser(0., 1.5*histMgenMC->GetMaximum());
  histMgenMC->SetLineColor(kRed+1);
  histMgenMC->SetLineColor(kRed+1);
  histMgenMC->SetMarkerColor(kRed+1);
  histMgenMC->SetMarkerStyle(5);
  histMgenMC->SetMarkerSize(0.75);
  histMgenMC->Draw("P X0 E1");
  
  divideByBinWidth(histMdetMC);
  histMdetMC->SetLineColor(kRed+1);
  histMdetMC->SetLineStyle(2);
  histMdetMC->Draw("SAME,HIST");
  
  divideByBinWidth(histMdetMCsig);
  histMdetMCsig->SetLineColor(kRed+1);
  histMdetMCsig->SetLineStyle(0);
  histMdetMCsig->Draw("SAME,HIST");
  
  divideByBinWidth(histMdetData);
  histMdetData->SetLineColor(kBlack);
  histMdetData->SetLineStyle(3);
  histMdetData->Draw("SAME,HIST");
  
  //divideByBinWidth(histMgenToy);
  //histMgenToy->SetLineColor(kMagenta+1);
  //histMgenToy->SetLineColor(kMagenta+1);
  //histMgenToy->SetMarkerColor(kMagenta+1);
  //histMgenToy->SetMarkerStyle(2);
  //histMgenToy->SetMarkerSize(0.75);
  //histMgenToy->Draw("SAME P X0 E1");
  
  divideByBinWidth(histMdetDataBGSubtracted);
  histMdetDataBGSubtracted->SetLineColor(kBlack);
  histMdetDataBGSubtracted->SetLineStyle(0);
  histMdetDataBGSubtracted->Draw("SAME,HIST");
  
  divideByBinWidth(histTotalError);
  histTotalError->SetLineColor(kBlack);
  histTotalError->SetMarkerColor(kBlack);
  histTotalError->SetMarkerStyle(20);
  histTotalError->SetMarkerSize(0.5);
  histTotalError->Draw("SAME P X0 E");
  divideByBinWidth(histMunfold);
  histMunfold->SetLineColor(kBlack);
  histMunfold->SetMarkerColor(kBlack);
  histMunfold->SetMarkerStyle(20);
  histMunfold->SetMarkerSize(0.5);
  histMunfold->Draw("SAME P X0 E1");
  
  divideByBinWidth(histMfoldback);
  histMfoldback->SetMarkerColor(kMagenta+1);
  histMfoldback->SetMarkerStyle(2);
  histMfoldback->SetMarkerSize(0.5);
  histMfoldback->Draw("SAME P X0 E1");
  
  TLegend leg(0.5,0.6,0.85,0.9);
  leg.SetLineWidth(0);
  leg.SetFillStyle(0);
  leg.AddEntry(histMgenMC, "MC gen", "p");
  leg.AddEntry(histMdetMC, "MC reco", "l");
  leg.AddEntry(histMdetMCsig, "MC reco signal", "l");
  //leg.AddEntry(histMgenToy, "toy gen", "p");
  leg.AddEntry(histMdetData, "data reco", "l");
  leg.AddEntry(histMdetDataBGSubtracted, "data (bg-sub)", "l");
  leg.AddEntry(histMunfold, "data unfolded", "ep");
  leg.AddEntry(histMfoldback, "data folded back", "p");
  leg.Draw();

  // show detector level distributions
  //    data (red)
  //    MC (black) [with completely wrong peak position and shape]
  //    unfolded data (blue)
  output.cd(3);
//  histMunfold->Draw("P X0 E1");
  //histMdetFold->SetLineColor(kBlue);
  //histMdetFold->Draw();
  //histMdetMC->Draw("SAME HIST");
  //
  //TH1 *histInput=unfold.GetInput("Minput",";mass(det)");

  //histInput->SetLineColor(kRed);
  //histInput->Draw("SAME");
  
  //TH2* ematrix = unfold.GetEmatrixInput("ematrix");
  gStyle->SetOptStat(0);
  histEmatTotal->Draw("colz,text");

  // show uncertainties
  output.cd(4);
  //histMdetGenToy->Draw("colz");
  //std::cout << histMdetGenToy->GetBinContent(5,5) << std::endl;
  
  std::vector<TH1F*> hists;
  bool first = true;
  int count = 0;
  TLegend legsys(0.,0.,1.,1.);
  for (TString uncertainty : uncertainties) {
    TH1* hist = unfold.GetDeltaSysSource(uncertainty, uncertainty);
    divideByBinWidth(hist);
    hist->Divide(histMunfold);
    hist->SetLineColor(++count);
    legsys.AddEntry(hist, uncertainty, "l");
    if (first) {
      hist->SetTitle("Systematic uncertainties");
      hist->GetYaxis()->SetRangeUser(-0.5,0.5);
      hist->Draw();
      first = false;
    }
    else hist->Draw("same");
  }
  for (int w = 1; w <= 20; ++w) {
    TH1* hist = unfold.GetDeltaSysSource(TString(std::to_string(w)), TString(std::to_string(w)));
    divideByBinWidth(hist);
    hist->Divide(histMunfold);
    hist->SetLineColor(++count);
    legsys.AddEntry(hist, TString("weight"+std::to_string(w)), "l");
    hist->Draw("same");
  }
  
  output.cd(5);
  legsys.Draw();
/*
  // show tau as a function of chi**2
  output.cd(5);
  //logTauX->Draw();
  rhoScan->Draw();
  rhoScan->SetTitle(";log_{10}(#tau);average(#rho_{i})");
  rhoScan->SetLineColor(kRed);
  bestRho->Draw("*");
  //bestLogTauLogChi2->SetMarkerColor(kRed);
  //bestLogTauLogChi2->Draw("*");
*/  

  for (int i = 1; i < histMunfold->GetNbinsX()+1; ++i) {
    std::cout << "dataUnfolded " << i << " " << histMunfold->GetBinContent(i) << " +/- " << histMunfold->GetBinError(i) << std::endl;
  }
  

  output.SaveAs("testUnfold0_"+observable+"_"+flavor+".pdf");
  output.SaveAs("testUnfold0_"+observable+"_"+flavor+".root");

  return 0;
}
