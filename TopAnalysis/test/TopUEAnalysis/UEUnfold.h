#ifndef UEUNFOLD_H
#define UEUNFOLD_H

#include <iostream>
#include <fstream>
#include <string.h>
#include <cstring>
#include <algorithm>
#include "Riostream.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TString.h"
#include "TUnfold.h"
#include "TUnfoldDensity.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <TH2D.h>
#include "TH1F.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaveText.h"
#include "TChain.h"
#include "TMath.h"
#include <TVector.h>
#include <TVector3.h>
#include "TCanvas.h"
#include "TFile.h"
#include "Riostream.h"
#include <fstream>
#include "TMatrixD.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TH2F.h"
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TF1.h>
#include <TStyle.h>
#include <TVector.h>
#include <TGraph.h>
#include "TUnfoldDensity.h"
#include "TSystem.h"

char name_h[100];
char file_name_2[100],name_h2[100], name_h3[100];
char file_axes[100],axes[2][100];
char name_toy_h[100];

bool toy_switch;

TH2F *m_incchmult;
TH2F *m_incchavgpt;
TH2F *m_incchflux;

TH1F *gen_incchmult;
TH1F *gen_incchavgpt;
TH1F *gen_incchflux;

TH1F *rec_incchmult;
TH1F *rec_incchavgpt;
TH1F *rec_incchflux;

TH1F *rec_incchmult_bkg[5];
TH1F *rec_incchavgpt_bkg[5];
TH1F *rec_incchflux_bkg[5];

TH1F *h_temp_incchmult[5];
TH1F *h_temp_incchavgpt[5];
TH1F *h_temp_incchflux[5];

TH1F *data_incchmult;
TH1F *data_incchavgpt;
TH1F *data_incchflux;


TH1F *rec_incchmult_toy;
TH1F *rec_incchavgpt_toy;
TH1F *rec_inccflux_toy;
 
TH1D *h_PULLS;
TH1D *h_PULLS_BINS[5];

TAxis *def_axis[2];

int nbins_gen_incchmult,nbins_rec_incchmult,nbinmax_gen_incchmult, nbinmax_rec_incchmult;
int cx_incchmult,cy_incchmult,cz_incchmult;

int nbins_gen_incchavgpt,nbins_rec_incchavgpt,nbinmax_gen_incchavgpt, nbinmax_rec_incchavgpt;
int cx_incchavgpt,cy_incchavgpt,cz_incchavgpt;

int nbins_gen_incchflux,nbins_rec_incchflux,nbinmax_gen_incchflux, nbinmax_rec_incchflux;
int cx_incchflux,cy_incchflux,cz_incchflux;


#endif
