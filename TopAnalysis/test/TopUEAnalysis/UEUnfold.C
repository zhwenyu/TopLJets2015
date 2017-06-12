#include "UEUnfold.h"

#include "TFile.h"
#include "TKey.h"
#include "TH2.h"
#include "TH1.h"
#include "TString.h"
#include "TStyle.h"

using namespace std;

TCanvas *dataControlCanvas(TH1 *totalData,TH1 *totaData_sub,TH1 *totalRec,TString tag="",TString extraLeg="");
TCanvas *showNormalizedMigrationMatrix(TH2 *migration);
TCanvas *showTauScan(TGraph *rhoScanGr,TGraph *bestRhoGr,double opt_tau,double opt_rho);
TObjArray *UEUnfold(char *dist,char *file_name,char *file_toy,char *sigName, float opt_tau,bool unfoldData);

//
TObjArray *UEUnfold(char *dist,char *file_name,char *file_toy,char *sigName,float opt_tau,bool unfoldData)
{
  TObjArray *results=new TObjArray();
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kTemperatureMap);
  gStyle->SetPaintTextFormat("4.0f");

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

  //migration matrix
  TH2 *totalMig   = (TH2 *)inF->Get( Form("%s_0_mig/%s_0_mig_%s",dist,dist,sigName) )->Clone("migration");
  totalMig->SetDirectory(0);
  inF->Close();

  results->Add( dataControlCanvas(totalData,totalDataSub,totalRec,"data") );
  results->Add( showNormalizedMigrationMatrix(totalMig) );

  //
  // READ TOY FILE
  //
  inF=TFile::Open(file_toy);

  //signal reconstructed
  TH1 *toyRec   = (TH1 *)inF->Get( Form("%s_None_True",dist) )->Clone("signal_toy");;
  toyRec->SetTitle("toy data");
  toyRec->SetDirectory(0);
  float sf=totalData->Integral()/toyRec->Integral();
  toyRec->Scale(sf);

  //subtract average contribution from fakes (use total fakes)
  TH1 *toyRecSub=(TH1 *)toyRec->Clone("signal_toy_sub");
  toyRecSub->SetTitle("toy data-fakes");
  toyRecSub->Add(totalFakes,-1);
  toyRecSub->SetDirectory(0);

  //MC truth for this toy
  TH1 *toyGen   = (TH1 *)inF->Get( Form("%s_None_False",dist) )->Clone("gen_toy");
  toyGen->SetTitle("toy MC truth");
  toyGen->SetDirectory(0);
  toyGen->Scale(sf); 
  inF->Close();
  
  results->Add( dataControlCanvas(toyRec,toyRecSub,totalRec,"toy") );

  //
  // UNFOLDING
  //

  //init TUnfoldDensity

  const double scaleBias(1.0);

  //TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
  TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
  
  //TUnfold::EConstraint constraint = TUnfold::kEConstraintArea;
  TUnfold::EConstraint constraint = TUnfold::kEConstraintNone;   
  
  TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;

  TUnfoldDensity unfold(totalMig,TUnfold::kHistMapOutputHoriz,regMode,constraint,densityFlags);  

  //
  // scan to determine optimal tau, if needed, and unfold the data
  //
  if(opt_tau<0)
    {       
      //unfold.SetInput(data,scaleBias);//test
      //   unfold.SetInput(totalRec,scaleBias);//WORKS EXACTLY
      unfold.SetInput(totalDataSub,scaleBias);

      Int_t nScan=100;
      double tauMin=1.e-4;
      double tauMax=1.e-1;
      TSpline *logTauX,*logTauY;
      TSpline *rhoScanSpline=0;

      //TUnfoldDensity::EScanTauMode tauflag = TUnfoldDensity::kEScanTauRhoAvg;
      //TUnfoldDensity::EScanTauMode tauflag = TUnfoldDensity::kEScanTauRhoMax;//crazy large correlations
      TUnfoldDensity::EScanTauMode tauflag = TUnfoldDensity::kEScanTauRhoSquareAvg;//crazy large correlations
      
      //unfold.ScanTau(nScan,tauMin,tauMax,&rhoScan,tauflag);	  
      int iBest = unfold.ScanTau(nScan,tauMin,tauMax,&rhoScanSpline,tauflag);	  

      // create graphs (rho scan and best tau)
      Double_t x,y;
      rhoScanSpline->GetKnot(iBest,x,y);
      opt_tau=unfold.GetTau();
      Double_t opt_rho=y;
      TGraph *bestRhoGr=new TGraph();
      bestRhoGr->SetPoint(0,x,y);
      TGraph *rhoScanGr=new TGraph();
      for(Int_t i=0;i<nScan;i++)
        {
          rhoScanSpline->GetKnot(i,x,y);
          rhoScanGr->SetPoint(i,x,y);
        }

      results->Add( showTauScan(rhoScanGr,bestRhoGr, opt_tau,opt_rho) );
      TVectorF *summary = new TVectorF(2);
      (*summary)[0]=opt_tau;
      (*summary)[1]=opt_rho;
      results->Add(summary);
    }

  //
  // unfold the data
  //
  if(unfoldData)
    {
      //unfold.DoUnfold(opt_tau,totalDataSub,scaleBias);

      TH1 *unfolded_data_raw = unfold.GetOutput("unfolded_data",0,0,"",kFALSE);
      TH1 *unfolded_data=(TH1 *)totalGen->Clone("corrected_data");
      unfolded_data->SetTitle("Data (corrected)");
      unfolded_data->SetDirectory(0);
      for(int xbin=1; xbin<=unfolded_data->GetNbinsX(); xbin++)
        {
          unfolded_data->SetBinContent(xbin,unfolded_data_raw->GetBinContent(xbin));
          unfolded_data->SetBinError(xbin,unfolded_data_raw->GetBinError(xbin));
        }

      TVectorF *summary = new TVectorF(3);
      (*summary)[0]=unfold.GetChi2A();
      (*summary)[1]=unfold.GetChi2L();
      (*summary)[2]=unfold.GetNdf();
      results->Add( summary );
      results->Add( unfolded_data );
      results->Add( dataControlCanvas(0,unfolded_data,totalGen,"unfolded",Form("#chi^{2}/ndf=(%3.1f+%3.1f)/%d",(*summary)[2],(*summary)[3],(int)(*summary)[4]) ) );
    }

  //
  // unfold toy
  //
  TUnfoldDensity toy_unfold(totalMig,TUnfold::kHistMapOutputHoriz,regMode,constraint,densityFlags);  
  toy_unfold.SetInput(toyRecSub,scaleBias);
  toy_unfold.DoUnfold(opt_tau);//,toyRecSub,scaleBias);
  TH1 *unfolded_toy_raw = toy_unfold.GetOutput("unfolded_signal_toy_sub",0,0,"",kFALSE);
  TH1 *unfolded_toy=(TH1 *)toyGen->Clone("corrected_toy");
  unfolded_toy->SetTitle("toy data (corrected)");
  unfolded_toy->SetDirectory(0);
  TH1 *bias=(TH1 *)toyGen->Clone("toy_bias");
  bias->SetTitle("bias");
  bias->Reset("ICE");
  TH1 *pull=(TH1 *)toyGen->Clone("toy_pull");
  pull->SetTitle("pull");
  pull->Reset("ICE");

  for(int xbin=1; xbin<=unfolded_toy_raw->GetNbinsX(); xbin++)
    {
      float unfVal(unfolded_toy_raw->GetBinContent(xbin)),unfUnc(unfolded_toy_raw->GetBinError(xbin));
      unfolded_toy->SetBinContent(xbin,unfVal);
      unfolded_toy->SetBinError(xbin,unfUnc);

      float targetVal(toyGen->GetBinContent(xbin));
      float delta(unfVal-targetVal);
      bias->SetBinContent(xbin,delta);
      if(unfUnc>0) pull->SetBinContent(xbin,delta/unfUnc);
    }
  results->Add( dataControlCanvas(0,unfolded_toy,toyGen,"toyunfolded") );
  results->Add(bias);
  results->Add(pull);

  return results;

  /*

   
   //data_bkg_sub->Scale(gen_incchmult->Integral()/data_bkg_sub->Integral());//gerekli diil gibi!
   unfold_incchmult.SetInput(data_incchmult_bkg_sub,scaleBias);
   unfold_incchavgpt.SetInput(data_incchavgpt_bkg_sub,scaleBias);
   unfold_incchflux.SetInput(data_incchflux_bkg_sub,scaleBias);

   unfold_incchmult.DoUnfold(opt_tau_incchmult);
   unfold_incchavgpt.DoUnfold(opt_tau_incchavgpt);
   unfold_incchflux.DoUnfold(opt_tau_incchflux);
   //unfold.DoUnfold(0);
   TH1 *unfolded_data_incchmult = unfold_incchmult.GetOutput("unfolded_data_incchmult",0,0,"",kFALSE);
   TH1 *unfolded_data_incchavgpt = unfold_incchavgpt.GetOutput("unfolded_data_incchavgpt",0,0,"",kFALSE);
   TH1 *unfolded_data_incchflux = unfold_incchflux.GetOutput("unfolded_data_incchflux",0,0,"",kFALSE);
   TH1 *folded_data_incchmult = unfold_incchmult.GetFoldedOutput("folded_data_incchmult",0,0,"",kFALSE);
   TH1 *folded_data_incchavgpt = unfold_incchavgpt.GetFoldedOutput("folded_data_incchavgpt",0,0,"",kFALSE);
   TH1 *folded_data_incchflux = unfold_incchflux.GetFoldedOutput("folded_data_incchflux",0,0,"",kFALSE);

   TH2 *EMatrixTotal_data_incchmult = unfold_incchmult.GetEmatrixTotal("unfolding total error matrix for incchmult");
   TH2 *EMatrixTotal_data_incchavgpt = unfold_incchavgpt.GetEmatrixTotal("unfolding total error matrix for incchavgpt");
   TH2 *EMatrixTotal_data_incchflux = unfold_incchflux.GetEmatrixTotal("unfolding total error matrix for incchflux");


   //unfold.SetInput(rec_incchmult_toy,scaleBias);

   // unfold.SetBias(gen_incchmult_real);

   //unfold.DoUnfold(0);//orig loc. 

   //unfold.SetInput(data_bkg_sub,scaleBias);//unfold data

   
   TCanvas *tcanv = new TCanvas();
   tcanv->Divide(2,2);
   tcanv->cd(1);
   rhoScan_incchmult->Draw();
   rhoScan_incchmult->SetTitle(";log_{10}(#tau);average(#rho_{i})");
   rhoScan_incchmult->SetLineColor(kRed);
   bestRho_incchmult->Draw("*");

   tcanv->cd(2);
   rhoScan_incchavgpt->Draw();
   rhoScan_incchavgpt->SetTitle(";log_{10}(#tau);average(#rho_{i})");
   rhoScan_incchavgpt->SetLineColor(kRed);
   bestRho_incchavgpt->Draw("*");

   tcanv->cd(3);
   rhoScan_incchflux->Draw();
   rhoScan_incchflux->SetTitle(";log_{10}(#tau);average(#rho_{i})");
   rhoScan_incchflux->SetLineColor(kRed);
   bestRho_incchflux->Draw("*");

   tcanv->SaveAs("knots.C");
   tcanv->SaveAs("knots.pdf");  
 
   TCanvas *c2 = new TCanvas();
   c2->Divide(2,2);
   //   float scale = rec_incchmult_real->GetSumOfWeights()/gen_incchmult_real->GetSumOfWeights();
   c2->cd(1);
   gen_incchmult->SetMarkerStyle(8);
   gen_incchmult->SetMarkerColor(2);
   gen_incchmult->SetLineColor(2);
   gen_incchmult->SetLineWidth(2); 
   gen_incchmult->Draw();

   // unfolded_incchmult->SetMarkerColor(2);
   // unfolded_incchmult->SetLineColor(2);
   // unfolded_incchmult->SetLineStyle(2);
   // unfolded_incchmult->SetMarkerStyle(8);
   // unfolded_incchmult->Draw("sames");

   unfolded_data_incchmult->SetMarkerStyle(21);
   unfolded_data_incchmult->SetMarkerColor(4);
   unfolded_data_incchmult->SetLineColor(4);
   unfolded_data_incchmult->SetLineWidth(2);
   unfolded_data_incchmult->Draw("sames");

   c2->cd(2);
   gen_incchavgpt->SetMarkerStyle(8);
   gen_incchavgpt->SetMarkerColor(2);
   gen_incchavgpt->SetLineColor(2);
   gen_incchavgpt->SetLineWidth(2); 
   gen_incchavgpt->Draw();

   // unfolded_incchavgpt->SetMarkerColor(2);
   // unfolded_incchavgpt->SetLineColor(2);
   // unfolded_incchavgpt->SetLineStyle(2);
   // unfolded_incchavgpt->SetMarkerStyle(8);
   // unfolded_incchavgpt->Draw("sames");

   unfolded_data_incchavgpt->SetMarkerStyle(21);
   unfolded_data_incchavgpt->SetMarkerColor(4);
   unfolded_data_incchavgpt->SetLineColor(4);
   unfolded_data_incchavgpt->SetLineWidth(2);
   unfolded_data_incchavgpt->Draw("sames");

   c2->cd(3);
   gen_incchflux->SetMarkerStyle(8);
   gen_incchflux->SetMarkerColor(2);
   gen_incchflux->SetLineColor(2);
   gen_incchflux->SetLineWidth(2); 
   gen_incchflux->Draw();

   // unfolded_incchflux->SetMarkerColor(2);
   // unfolded_incchflux->SetLineColor(2);
   // unfolded_incchflux->SetLineStyle(2);
   // unfolded_incchflux->SetMarkerStyle(8);
   // unfolded_incchflux->Draw("sames");

   unfolded_data_incchflux->SetMarkerStyle(21);
   unfolded_data_incchflux->SetMarkerColor(4);
   unfolded_data_incchflux->SetLineColor(4);
   unfolded_data_incchflux->SetLineWidth(2);
   unfolded_data_incchflux->Draw("sames");

   c2->Update();
   c2->SaveAs("unfolded.C");
   
   cout<<"before pulls"<<endl;
   ofstream pulls_file,pulls_file_bin[5];
   pulls_file.open("pulls_file.txt",ios::app);

   char namepullbinfile[100];
   for (int i =0;i<5;i++){ 
     sprintf(namepullbinfile,"pulls_file_bin%i.txt",i);
     pulls_file_bin[i].open(namepullbinfile,ios::app);
   }
   h_PULLS = new TH1D("h_PULLS","h_PULLS",100,-5,5);
  
   float pull = -99.;
   float pull_b[5];
   TCanvas *c_pull = new TCanvas();
   c_pull->cd();
   for (int i=0;i<5;i++){ 
     //if (unfolded->GetBinError(i+1)) 
     pull_b[i] = -999.;
     pull = (unfolded_toy_incchmult->GetBinContent(i+1)-gen_incchmult->GetBinContent(i+1))/unfolded_toy_incchmult->GetBinError(i+1);
     pull_b[i] = (unfolded_toy_incchmult->GetBinContent(i+1)-gen_incchmult->GetBinContent(i+1))/unfolded_toy_incchmult->GetBinError(i+1);
     h_PULLS->Fill(pull);
     cout<<"pull toy, gen toy, err:  "<<unfolded_toy_incchmult->GetBinContent(i+1)<<"  "<<gen_incchmult->GetBinContent(i+1)<<"  "<<unfolded_toy_incchmult->GetBinError(i+1)<<endl;
     cout<<"pull:   "<<pull<<endl;
     pulls_file<< pull <<"\n";
     pulls_file_bin[i]<<pull_b[i]<<"\n";
   }
   pulls_file.close();
   for (int i =0;i<5;i++) pulls_file_bin[i].close();
   h_PULLS->Draw();
   c_pull->SaveAs("pull.C");

  
   TCanvas *c4 = new TCanvas();
   c4->Divide(2,2);
   c4->cd(1);
   data_incchmult_bkg_sub->Draw();
   folded_data_incchmult->SetLineColor(4);
   folded_data_incchmult->SetMarkerColor(4);
   folded_data_incchmult->Draw("sames");

   c4->cd(2);
   data_incchavgpt_bkg_sub->Draw();
   folded_data_incchavgpt->SetLineColor(4);
   folded_data_incchavgpt->SetMarkerColor(4);
   folded_data_incchavgpt->Draw("sames");

   c4->cd(3);
   data_incchflux_bkg_sub->Draw();
   folded_data_incchflux->SetLineColor(4);
   folded_data_incchflux->SetMarkerColor(4);
   folded_data_incchflux->Draw("sames");
   c4->SaveAs("folded.C");

  TCanvas *c5 = new TCanvas();
  c5->cd(1);
  gen_incchmult->SetLineColor(1);
  gen_incchmult->SetLineWidth(2);
  gen_incchmult->Draw("e1");
  unfolded_incchmult->SetMarkerColor(2);
  unfolded_incchmult->SetLineColor(2);
  unfolded_incchmult->SetLineStyle(2);
  unfolded_incchmult->SetMarkerStyle(8);
  unfolded_incchmult->Draw("sames");
  unfolded_data_incchmult->SetMarkerColor(4);
  unfolded_data_incchmult->SetLineColor(4);
  unfolded_data_incchmult->SetLineStyle(3);
  unfolded_data_incchmult->Draw("sames");
  c5->SaveAs("unfolded_incmult_and_sqrtN.C");

  TCanvas *c6 = new TCanvas();
  c6->Divide(2,2);
  c6->cd(1);
  EMatrixTotal_data_incchmult->Draw("colz");
  c6->cd(2);
  EMatrixTotal_data_incchavgpt->Draw("colz");
  c6->cd(3);
  EMatrixTotal_data_incchflux->Draw("colz");
  c6->SaveAs("error_matrices.C");
  c6->SaveAs("error_matrices.pdf");

*/
}

TCanvas *dataControlCanvas(TH1 *totalData,TH1 *totalDataSub,TH1 *totalRec,TString tag,TString extraLeg)
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
TCanvas *showNormalizedMigrationMatrix(TH2 *migration)
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
TCanvas *showTauScan(TGraph *rhoScanGr,TGraph *bestRhoGr,double opt_tau,double opt_rho)
{
  TCanvas *c = new TCanvas("cscan","cscan",500,500);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetBottomMargin(0.1);
  c->SetLeftMargin(0.1);

  rhoScanGr->Draw("ac");
  rhoScanGr->GetXaxis()->SetTitle("log_{10}(#tau)");
  rhoScanGr->GetYaxis()->SetTitle("Global #rho");
  bestRhoGr->SetMarkerStyle(29);
  bestRhoGr->SetMarkerColor(kRed);
  bestRhoGr->SetLineColor(kRed);
  bestRhoGr->SetMarkerSize(1.5);
  bestRhoGr->Draw("p");

  TLatex *txt=new TLatex();
  txt->SetTextFont(42);
  txt->SetTextSize(0.04);
  txt->SetNDC();
  txt->DrawLatex(0.12,0.96,"#bf{CMS} #it{Simulation Preliminary}");
  txt->DrawLatex(0.7,0.97,"#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
  txt->DrawLatex(0.12,0.9,Form("#tau=%3.4f #rho=%3.3f",opt_tau,opt_rho));
  return c;

}
