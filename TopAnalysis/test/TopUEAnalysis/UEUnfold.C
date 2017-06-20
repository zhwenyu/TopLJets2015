#include "UEUnfold.h"

#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TH2.h"
#include "TH1.h"
#include "TString.h"
#include "TStyle.h"

using namespace std;

class UEUnfold
{
public:
  UEUnfold() : results_(0)
  { 
    initROOTStyle(); 
    reset(); 
  }
  void unfoldData(char *dist,char *file_name,char *syst_file_name,char *sigName);
  void unfoldToy(char *dist,char *file_toy,float opt_tau, TH2 *mig,float norm,TH1 *fakes);
  float doUnfold(float opt_tau, TH2 *mig, TH1 *data,TH1 *gen,bool storeFull=true,TString pfix="");
  TObjArray *getResults() { return results_; }
  void reset() 
  { 
    if(results_) results_->Delete(); 
    results_=new TObjArray; 
  }
  ~UEUnfold() { };

private:
  void initROOTStyle();
  float performTauScan(TH2 *mig, TH1 *data,TUnfoldDensity &unfold);

  TCanvas *ratioCanvas(TH1 *num,TH1 *den,TString name, TString xtitle,TString ytitle,Float_t ymin,Float_t ymax);
  TCanvas *dataControlCanvas(TH1 *totalData,TH1 *totaData_sub,TH1 *totalRec,TString tag="",TString extraLeg="");
  TCanvas *showNormalizedMigrationMatrix(TH2 *migration);
  TCanvas *showTauScan(TGraph *rhoScanGr,TGraph *bestRhoGr, TGraph *rho2ScanGr,double opt_tau,double opt_rho);
  TH1 * bottomLineTest(TH1 *totalDataSub,TH1 *toyRecSub, TH1*unfolded_data, TH1 *toyGen, TH2 *unfolded_ematTotal,bool doNorm=true);
  TObjArray *results_;
};


//
void UEUnfold::unfoldData(char *dist,char *file_name,char *syst_file_name,char *sigName)
{
  //
  // READ FILE WITH SUMMED UP CONTRIBUTIONS
  //
  TFile *inF=TFile::Open(file_name);

  //reconstruction level distributions
  TH1 *totalRec=0,*totalData=0,*totalBkg=0;
  TString recDir(Form("%s_None_True",dist));
  TIter next( ((TDirectory *)inF->Get(recDir))->GetListOfKeys() );
  TKey *key;
  while ((key = (TKey*)next())) 
    {
      TObject *obj=key->ReadObj();
      if (!obj->InheritsFrom("TH1")) continue;
      TH1 *h = (TH1*)obj;
      TString title(h->GetTitle());
      if(title=="Data")       { totalData=(TH1 *)h->Clone("data"); totalData->SetDirectory(0); }
      else if(title==sigName) { totalRec=(TH1 *)h->Clone("signal"); totalRec->SetDirectory(0); }
      else {
        if(totalBkg==0) { totalBkg=(TH1 *)h->Clone("bkg"); totalBkg->SetDirectory(0); }
        else            { totalBkg->Add(h); }
      }
    }

  totalRec->SetTitle("signal");
  totalBkg->SetTitle("bkg");

  //fake contributions to the data
  TH1 *totalFakes = (TH1 *)inF->Get( Form("%s_fakes_True/%s_fakes_True_%s",dist,dist,sigName) )->Clone("fakes");
  totalFakes->SetTitle("fakes");
  totalFakes->SetDirectory(0);

  //subtract backgrounds and fakes from data
  TH1 *totalDataSub = (TH1 *)totalData->Clone("datasub");
  totalDataSub->SetTitle("Data-bkg-fakes");
  totalDataSub->Add(totalFakes,-1);
  totalDataSub->Add(totalBkg,-1);
  totalDataSub->SetDirectory(0);

  //gen level distribution
  TH1 *totalGen   = (TH1 *)inF->Get( Form("%s_None_False/%s_None_False_%s",dist,dist,sigName) )->Clone("gen");
  totalGen->SetTitle("Gen. level");
  totalGen->SetDirectory(0);

  //migration matrices (0 is the nominal)
  std::map<TString, TH2 *> migMatrices;
  for(int i=0; i<=25; i++)
    {
      TString pfix( i==0 ? "" : Form("_%d",i) );
      migMatrices[pfix]=(TH2 *)inF->Get( Form("%s_%d_mig/%s_%d_mig_%s",dist,i,dist,i,sigName) )->Clone("migration"+pfix);
      migMatrices[pfix]->SetDirectory(0);
    }
  inF->Close();
  
  //add simulated systematics
  TFile *systF=TFile::Open(syst_file_name);
  TList *keys=((TDirectory *)systF->Get(Form("%s_0_mig",dist)))->GetListOfKeys();
  TIter nextKey(keys);
  while ((key = (TKey*)nextKey()))
    {
      TH2 *h=(TH2 *)key->ReadObj();
      if(h->Integral()==0) continue;
      TString name=h->GetName();
      if(!name.Contains(sigName) || name.Contains("t#bar{t}t#bar{t}")) continue;
      TString pfix=h->GetTitle(); 
      migMatrices[pfix]=(TH2 *)h->Clone("migration"+pfix);
      migMatrices[pfix]->Scale(migMatrices[pfix]->Integral()/migMatrices[""]->Integral());
      migMatrices[pfix]->SetDirectory(0);
    }
  systF->Close();

  //save results
  results_->Add(totalDataSub);
  results_->Add(migMatrices[""]);
  results_->Add(totalFakes);
  results_->Add( ratioCanvas(totalDataSub,totalData,"datasub","Reconstructed bin","(Data-Fakes-Background)/Data",0.8,1.0) );
  results_->Add( dataControlCanvas(totalData,totalDataSub,totalRec,"data") );
  results_->Add( showNormalizedMigrationMatrix(migMatrices[""]) );
  float opt_tau=doUnfold(-1,migMatrices[""],totalDataSub,totalGen);
  cout << "After unfolding with nominal matrix tau=" << opt_tau << endl;
  for(std::map<TString, TH2 *>::iterator it=migMatrices.begin(); it!=migMatrices.end(); it++)
    {
      if(it->first=="") continue;
      cout << "Unfolding variation" << it->first << endl;
      doUnfold(opt_tau,it->second,totalDataSub,totalGen,false,it->first);
    }
  doUnfold(opt_tau,migMatrices[""],totalData,totalGen,false,"bckpfakes");
}

//
float UEUnfold::doUnfold(float opt_tau, TH2 *mig, TH1 *data,TH1 *gen,bool storeFull,TString pfix)
{
  const double scaleBias(1.0);

  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;                            //regularization mode
  //TUnfold::ERegMode regMode = TUnfold::kRegModeSize;

  TUnfold::EConstraint constraint = TUnfold::kEConstraintNone;                       //area constraint  
  //TUnfold::EConstraint constraint = TUnfold::kEConstraintArea;                    

  TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeeNone;       //bin definition  
  //TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;  

  TUnfoldDensity unfold(mig,TUnfold::kHistMapOutputHoriz,regMode,constraint,densityFlags);  
  unfold.SetInput(data,scaleBias);
 
  // scan to determine optimal tau, if needed
  if(opt_tau<0) opt_tau=performTauScan(mig,data,unfold);

  //unfold with optimized tau  
  unfold.DoUnfold(opt_tau);

  //get results
  TH1 *folded_unfolded_data = unfold.GetFoldedOutput("folded_data",0,0,"",kFALSE);
  TH1 *unfolded_data_raw = unfold.GetOutput("unfolded_data",0,0,"",kFALSE);
  TH1 *unfolded_data=(TH1 *) gen->Clone("corrected_data"+pfix);
  unfolded_data->SetTitle("Data (corrected)");
  unfolded_data->SetDirectory(0);
  TH2 *unfolded_ematTotal = unfold.GetEmatrixTotal("EmatTotal");
  for(int xbin=1; xbin<=unfolded_data->GetNbinsX(); xbin++)
    {
      unfolded_data->SetBinContent(xbin,unfolded_data_raw->GetBinContent(xbin));
      unfolded_data->SetBinError(xbin,TMath::Sqrt(unfolded_ematTotal->GetBinContent(xbin,xbin)));
    }

  //save summary
  if(storeFull)
    {
      TVectorF *summary = new TVectorF(3);
      (*summary)[0]=unfold.GetChi2A();
      (*summary)[1]=unfold.GetChi2L();
      (*summary)[2]=unfold.GetNdf();
      results_->Add( summary );
      results_->Add( unfolded_ematTotal );
      results_->Add( folded_unfolded_data->Clone("data_particle_folded") );
      results_->Add( dataControlCanvas(0,unfolded_data,gen,"unfolded",Form("#chi^{2}/ndf=(%3.1f+%3.1f)/%d",(*summary)[0],(*summary)[1],(int)(*summary)[2]) ) );
      results_->Add( ratioCanvas(folded_unfolded_data,data,"datafold","Reconstructed bin","Folded/Reconstructed",0.8,1.2) );
      results_->Add( unfolded_data_raw );
    }
  results_->Add( unfolded_data );

  return opt_tau;
} 
 




//
float UEUnfold::performTauScan(TH2 *mig, TH1 *data,TUnfoldDensity &unfold)
{  
  Int_t nScan=100;
  double tauMin=1.e-4;
  double tauMax=1.e-1;
  TSpline *logTauX,*logTauY;
  TSpline *rhoScanSpline=0;
  
  TUnfoldDensity::EScanTauMode tauflag = TUnfoldDensity::kEScanTauRhoAvg;          //scan mode
  //TUnfoldDensity::EScanTauMode tauflag = TUnfoldDensity::kEScanTauRhoMax;
  //TUnfoldDensity::EScanTauMode tauflag = TUnfoldDensity::kEScanTauRhoSquareAvg;
  
  int iBest = unfold.ScanTau(nScan,tauMin,tauMax,&rhoScanSpline,tauflag);	  
  
  // create graphs (rho scan and best tau)
  Double_t x,y;
  rhoScanSpline->GetKnot(iBest,x,y);
  float opt_tau=unfold.GetTau();
  Double_t opt_rho=y;
  TGraph *bestRhoGr=new TGraph();
  bestRhoGr->SetPoint(0,x,y);
  TGraph *rhoScanGr=new TGraph();
  for(Int_t i=0;i<nScan;i++)
    {
      rhoScanSpline->GetKnot(i,x,y);
      rhoScanGr->SetPoint(i,x,y);
    }
  
  //for reference display also the correlation squared based optimization
  unfold.ScanTau(nScan,tauMin,tauMax,&rhoScanSpline,TUnfoldDensity::kEScanTauRhoSquareAvg);
  TGraph *rho2ScanGr=new TGraph();
  for(Int_t i=0;i<nScan;i++)
    {
      rhoScanSpline->GetKnot(i,x,y);
      rho2ScanGr->SetPoint(i,x,sqrt(y));
    }
  
  //save results
  results_->Add( showTauScan(rhoScanGr,bestRhoGr, rho2ScanGr, opt_tau,opt_rho) );  
  TVectorF *summary = new TVectorF(2);
  (*summary)[0]=opt_tau;
  (*summary)[1]=opt_rho;
  results_->Add(summary);

  return opt_tau;
}

                    

//
void UEUnfold::unfoldToy(char *dist,char *file_toy,float opt_tau,TH2 *mig,float norm,TH1 *fakes)
{
  //
  // READ TOY FILE
  //
  TFile *inF=TFile::Open(file_toy);

  //signal reconstructed
  TH1 *toyRec   = (TH1 *)inF->Get( Form("%s_None_True",dist) )->Clone("signal_toy");;
  toyRec->SetTitle("toy data");
  toyRec->SetDirectory(0);
  float sf=norm/toyRec->Integral();
  toyRec->Scale(sf);

  //subtract average contribution from fakes (use total fakes)
  TH1 *toyRecSub=(TH1 *)toyRec->Clone("signal_toy_sub");
  toyRecSub->SetTitle("toy data-fakes");
  if(fakes) toyRecSub->Add(fakes,-1);
  toyRecSub->SetDirectory(0);

  //MC truth for this toy
  TH1 *toyGen   = (TH1 *)inF->Get( Form("%s_None_False",dist) )->Clone("gen_toy");
  toyGen->SetTitle("toy MC truth");
  toyGen->SetDirectory(0);
  toyGen->Scale(sf); 
  inF->Close();
  
  doUnfold(opt_tau,mig,toyRecSub,toyGen);

  TH1 *unfolded_toy=0;
  for(int i=0; i<results_->GetEntriesFast(); i++)
    {
      TString rname=results_->At(i)->GetName();
      if(rname!="corrected_data") continue;
      unfolded_toy=(TH1 *)results_->At(i);
    }

  //
  // unfold toy
  //
  TH1 *bias=(TH1 *)toyGen->Clone("toy_bias");
  bias->SetTitle("bias");
  bias->Reset("ICE");
  TH1 *pull=(TH1 *)toyGen->Clone("toy_pull");
  pull->SetTitle("pull");
  pull->Reset("ICE");
  for(int xbin=1; xbin<=unfolded_toy->GetNbinsX(); xbin++)
    {
      float unfVal(unfolded_toy->GetBinContent(xbin));
      float unfUnc(unfolded_toy->GetBinError(xbin));
      float targetVal(toyGen->GetBinContent(xbin));
      float targetUnc(toyGen->GetBinError(xbin));
      float delta(unfVal-targetVal);
      bias->SetBinContent(xbin,delta);
      if(unfUnc>0) pull->SetBinContent(xbin,delta/unfUnc);
    }
  results_->Add( dataControlCanvas(0,unfolded_toy,toyGen,"toyunfolded") );
  results_->Add(bias);
  results_->Add(pull);
  
  //  results_->Add(bottomLineTest(totalDataSub,toyRecSub,unfolded_data,toyGen,unfolded_ematTotal));

}

//
TH1 *UEUnfold::bottomLineTest(TH1 *totalDataSub,TH1 *toyRecSub, TH1*unfolded_data, TH1 *toyGen, TH2 *unfolded_ematTotal,bool doNorm)
{
  /*
  //
  // bottom line test
  //
  
  //smeared level
  int nbinsx(totalDataSub->GetNbinsX());
  TVectorF smearedDiff(nbinsx);
  TMatrixF smearedCov(nbinsx,nbinsx);
  float normData(totalDataSub->Integral()), normToy(toyRecSub->Integral());
  for(int xbin=1; xbin<=nbinsx; xbin++)
    {
      if(doNorm)
        {
          smearedDiff[xbin-1]=totalDataSub->GetBinContent(xbin)/normData-toyRecSub->GetBinContent(xbin)/normToy;
          smearedCov[xbin-1][xbin-1]=pow(totalDataSub->GetBinError(xbin)/normData,2)+pow(toyRecSub->GetBinError(xbin)/normToy,2);
        }
      else
        {
          smearedDiff[xbin-1]=totalDataSub->GetBinContent(xbin)-toyRecSub->GetBinContent(xbin);
          smearedCov[xbin-1][xbin-1]=pow(totalDataSub->GetBinError(xbin),2)+pow(toyRecSub->GetBinError(xbin),2);
        }
    }
  TVectorF a=smearedDiff;
  a *= smearedCov.Invert();
  float smearedChi2=a*smearedDiff;
  float pval=TMath::Prob(smearedChi2,nbinsx);
  
  //unfolded level
  int nbinsxpart=unfolded_data->GetNbinsX();
  TVectorF particleDiff(nbinsxpart);
  TMatrixF particleCov(nbinsxpart,nbinsxpart);
  float normUnfData(unfolded_data->Integral()), normToyGen(toyGen->Integral());
  for(int xbin=1; xbin<=nbinsxpart; xbin++)
    {
      if(doNorm)
        {
          particleDiff[xbin-1]=unfolded_data->GetBinContent(xbin)/normUnfData-toyGen->GetBinContent(xbin)/normToyGen;
          for(int ybin=1; ybin<=unfolded_ematTotal->GetNbinsY(); ybin++)
            {
              particleCov[xbin-1][ybin-1]=unfolded_ematTotal->GetBinContent(xbin,ybin)/pow(normUnfData,2);
              if(xbin==ybin) particleCov[xbin-1][ybin-1] += pow(toyGen->GetBinError(xbin)/normToyGen,2);
            }
        }
      else
        {
          particleDiff[xbin-1]=unfolded_data->GetBinContent(xbin)-toyGen->GetBinContent(xbin);
          for(int ybin=1; ybin<=unfolded_ematTotal->GetNbinsY(); ybin++)
            {
              particleCov[xbin-1][ybin-1]=unfolded_ematTotal->GetBinContent(xbin,ybin);
              if(xbin==ybin) particleCov[xbin-1][ybin-1] += pow(toyGen->GetBinError(xbin),2);
            }
        }
    }
  TVector b=particleDiff;
  b *= particleCov.Invert();
  float particleChi2=b*particleDiff;
  float particlepval=TMath::Prob(particleChi2,nbinsxpart);
  
  cout << "Smeared\t\tUnfolded" << endl;
  cout << smearedChi2 << " " << pval << " | " << particleChi2 << " " << particlepval << endl;
  cout << toyRecSub->Chi2Test(totalDataSub) << " " << toyGen->Chi2Test(unfolded_data) << endl;
  */

  TH1 *bottomLine=new TH1F("bottomline","",2,0,2);
  //  bottomLine->SetBinContent(1,smearedChi2/nbinsx);
  //bottomLine->SetBinContent(2,pval);
  //bottomLine->SetBinContent(3,particleChi2/nbinsxpart);
  // bottomLine->SetBinContent(4,particlepval);
  bottomLine->SetDirectory(0);
  return bottomLine;
}

//
TCanvas *UEUnfold::ratioCanvas(TH1 *num,TH1 *den,TString name, TString xtitle,TString ytitle,Float_t ymin,Float_t ymax)
{
  TCanvas *c = new TCanvas("c"+name,"c"+name,500,500);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.02);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.1);
  c->SetGridy();
  TH1 *ratio=(TH1 *)num->Clone("ratio");
  ratio->Divide(den);
  ratio->SetMarkerStyle(20);
  ratio->Draw("ep");
  ratio->GetYaxis()->SetTitle(ytitle);
  ratio->GetXaxis()->SetTitle(xtitle);
  ratio->GetYaxis()->SetRangeUser(ymin,ymax); 
  TLatex *txt=new TLatex();
  txt->SetTextFont(42);
  txt->SetTextSize(0.04);
  txt->SetNDC();
  txt->DrawLatex(0.15,0.9,"#bf{CMS} #it{Preliminary}");
  txt->DrawLatex(0.75,0.97,"#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
  return c;

}

TCanvas *UEUnfold::dataControlCanvas(TH1 *totalData,TH1 *totalDataSub,TH1 *totalRec,TString tag,TString extraLeg)
{
  TCanvas *c = new TCanvas("c"+tag,"c"+tag,500,500);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.02);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.1);
  if(totalData) { totalData->SetMarkerStyle(20);    totalData->SetMarkerColor(kGray); totalData->SetLineColor(kGray); }
  totalDataSub->SetMarkerStyle(20); totalDataSub->SetMarkerColor(1);  totalDataSub->SetLineColor(1);
  totalRec->SetFillStyle(0);        totalRec->SetLineColor(1); 
  totalDataSub->Draw("p");
  if(totalData) totalData->Draw("psame");
  totalRec->Draw("histsame");
  TLegend *leg=new TLegend(0.15,0.9,0.4,0.7);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->AddEntry(totalDataSub,totalDataSub->GetTitle(),"ep");
  if(totalData) leg->AddEntry(totalData,totalData->GetTitle(),"ep");
  leg->AddEntry(totalRec,totalRec->GetTitle(),"l");
  leg->Draw();
  TLatex *txt=new TLatex();
  txt->SetTextFont(42);
  txt->SetTextSize(0.04);
  txt->SetNDC();
  txt->DrawLatex(0.15,0.9,"#bf{CMS} #it{Preliminary}");
  txt->DrawLatex(0.5,0.9,"#scale[0.8]{#it{"+tag+"}}");
  if(extraLeg!="")   txt->DrawLatex(0.5,0.85,extraLeg);
  txt->DrawLatex(0.75,0.97,"#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
  return c;
}

//
TCanvas *UEUnfold::showNormalizedMigrationMatrix(TH2 *migration)
{
  //normalize to 100% by columns
  TH2* normmig=(TH2 *)migration->Clone("migration_norm");
  normmig->SetDirectory(0);
  for(int xbin=1; xbin<=migration->GetNbinsX(); xbin++)
    {
      TH1 *tmp=migration->ProjectionY("tmp",xbin,xbin);
      float total=tmp->Integral(1,tmp->GetNbinsX());
      tmp->Delete();
      if(total==0) continue;
      for(int ybin=1; ybin<=migration->GetNbinsY(); ybin++)
        {
          float val=migration->GetBinContent(xbin,ybin);
          float unc=migration->GetBinError(xbin,ybin);
          normmig->SetBinContent(xbin,ybin,100.*val/total);
          normmig->SetBinError(xbin,ybin,100.*unc/total);
        }
    }
  
  //show it
  TCanvas *c = new TCanvas("cmig","cmig",500,500);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.15);
  c->SetBottomMargin(0.1);
  c->SetLeftMargin(0.1);
  normmig->Draw("colz text");
  normmig->GetZaxis()->SetLabelSize(0.03);
  normmig->GetYaxis()->SetLabelSize(0.03);
  normmig->GetXaxis()->SetLabelSize(0.03);
  normmig->GetXaxis()->SetTitle("Generator level bin");
  normmig->GetYaxis()->SetTitle("Reconstruction level bin");
  normmig->GetZaxis()->SetRangeUser(0.,100.);
  TLatex *txt=new TLatex();
  txt->SetTextFont(42);
  txt->SetTextSize(0.04);
  txt->SetNDC();
  txt->DrawLatex(0.12,0.96,"#bf{CMS} #it{Simulation Preliminary}");
  txt->DrawLatex(0.7,0.97,"#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
  return c;
 
}

//
TCanvas *UEUnfold::showTauScan(TGraph *rhoScanGr,TGraph *bestRhoGr,TGraph *rho2ScanGr,double opt_tau,double opt_rho)
{
  TCanvas *c = new TCanvas("cscan","cscan",500,500);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetBottomMargin(0.1);
  c->SetLeftMargin(0.1);

  rhoScanGr->Draw("ac");
  rhoScanGr->GetXaxis()->SetTitle("log_{10}(#tau)");
  rhoScanGr->GetYaxis()->SetTitle("<#rho> or #sqrt{<#rho^{2}>}");
  rho2ScanGr->Draw("c");
  rho2ScanGr->SetLineColor(kGray);
  bestRhoGr->SetMarkerStyle(29);
  bestRhoGr->SetMarkerColor(kRed);
  bestRhoGr->SetLineColor(kRed);
  bestRhoGr->SetMarkerSize(1.5);
  bestRhoGr->Draw("p");

  TLatex *txt=new TLatex();
  txt->SetTextFont(42);
  txt->SetTextSize(0.04);
  txt->SetNDC();
  txt->DrawLatex(0.12,0.96,"#bf{CMS} #it{preliminary}");
  txt->DrawLatex(0.7,0.97,"#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
  txt->DrawLatex(0.12,0.9,Form("#tau=%3.4f #rho=%3.3f",opt_tau,opt_rho));
  return c;

}

//
void UEUnfold::initROOTStyle()
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kTemperatureMap);
  gStyle->SetPaintTextFormat("4.0f");
}
