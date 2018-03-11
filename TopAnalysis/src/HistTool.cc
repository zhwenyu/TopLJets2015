#include "TopLJets2015/TopAnalysis/interface/HistTool.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TH2F.h"

// Histogram tool for automatic creation of 2D uncertainty histograms
HistTool::HistTool(unsigned int nsyst) :
  nsyst_(nsyst)
{}

//
void HistTool::addHist(TString title, TH1* hist) {
  if(hist->InheritsFrom("TH2")) {
    all2dPlots_[title]=(TH2 *)hist;
  }
  else {
    allPlots_[title] = hist;
    if (nsyst_ > 0) {
      all2dPlots_[title+"_syst"] = new TH2F(title+"_syst", hist->GetTitle(), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), nsyst_+1, -0.5, nsyst_+0.5);
      all2dPlots_[title+"_syst"]->SetXTitle(hist->GetXaxis()->GetTitle());
      all2dPlots_[title+"_syst"]->SetYTitle("Variation (0=default)");
    }
  }
}

//
void HistTool::fill(TString title, double value, std::vector<double> weights) {
  if (not allPlots_.count(title)) {
    std::cout << "Histogram " << title << " not registered, not filling." << std::endl;
    return;
  }
  if(allPlots_[title]->InheritsFrom("TH2")) return;
  allPlots_[title]->Fill(value, weights[0]);
  if (nsyst_ > 0) {
    if (weights.size() > nsyst_)
      std::cout << "WARNING: Size of uncertainty weight vector larger than uncertainty histogram size." << std::endl;
    all2dPlots_[title+"_syst"]->Fill(value, 0., weights[0]);
    for (unsigned int i = 1; i < weights.size(); ++i) {
      all2dPlots_[title+"_syst"]->Fill(value, i, weights[0]*weights[i]);
    }
  }
}


