#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/MethodCategory.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"


//
void RunBjetChargeClassification(TString inname,
                                 TString cut="nch>0 && vtxnch==0 && (nmu+vtxnmu)==0",
                                 TString MethodName="BDT_Simple")
{
   // This loads the library
   TMVA::Tools::Instance();

   //prepare output
   TString outname(MethodName+".root");
   TFile* outputFile = TFile::Open( outname, "RECREATE" );

   //init TMVA
   std::string factoryOptions( "!V:!Silent:Transformations=I;D;P;G,D:!Color:!DrawProgressBar" );
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassificationCategory", outputFile, factoryOptions );

   // Define the input variables used for the MVA training
   factory->AddVariable( "abseta := fabs(eta)",     'F' );
   if(cut.Contains("nch>0"))
     {
       factory->AddVariable( "nch",      'I' );
       factory->AddVariable( "ch",       'F' );
       factory->AddVariable( "nmu",      'I' );
       factory->AddVariable( "much",     'F' );
     }
   if(cut.Contains("vtxnch>0"))
     {
       factory->AddVariable( "vtxnch",   'I' );
       factory->AddVariable( "vtxch",    'F' );
       factory->AddVariable( "vtxnmu",   'I' );
       factory->AddVariable( "vtxmuch",  'F' );
       factory->AddVariable( "vtxmass",    'F' );
       factory->AddVariable( "vtxchi2",    'F' );
       factory->AddVariable( "vtxL3d",    'F' );
       factory->AddVariable( "vtxL3dSig",    'F' );

     }
   
   // Load the signal and background event samples from ROOT trees
   TFile *input = TFile::Open( inname );
   
   //add the data for the training
   TCut mycuts("g_bId>0 &&" + cut);
   TCut mycutb("g_bId<0 &&" +cut);
   TTree *data     = (TTree*)input->Get("data");
   factory->AddSignalTree( data,     1.0);
   factory->AddBackgroundTree( data, 1.0);
   factory->PrepareTrainingAndTestTree( mycuts, 
                                        mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // --- Categorised classifier
   TString cfg(
               "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"
               );
   factory->BookMethod( TMVA::Types::kBDT,MethodName,cfg);


   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // Save the output
   outputFile->Close();

   // Clean up
   delete factory;
}

//
void BjetChargeClassification(TString inname)
{
  TString cuts[]={"nch>0 && vtxnch==0",
                  "nch>0 && vtxnch>0"};

  TString methods[]={"BDT_Simple",
                     "BDT_Vtx"};

  for(int i=0; i<2; i++)
    RunBjetChargeClassification(inname,cuts[i],methods[i]);
}
