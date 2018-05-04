#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/MethodBDT.h"

#include "TMVA/ResultsClassification.h"
#include "TMVA/ROCCurve.h"
// Authored by Hamed Bakhshiansohi 08/12/2016
// I edited tmva/tmvagui/src/mvas.cxx file
// and changed "x" -> "d" option for kolmogorov test so that the overtraining works
// you need to compile the root again

using namespace std;
using namespace TMVA;

class BDTOptions{
public:
  int ntrees ;
  double minnodesize ;
  int maxdepth ;
  double adaboostbeta ;
  int ncuts ;
  TString name;
  
  BDTOptions(  TString name_ = "" , int ntrees_ = 50 , double minnodesize_ = 2. , int maxdepth_ = 50 , double adaboostbeta_ = 0.5 , int ncuts_ = 10  );
  void Set( string options_ );
  TString to_string( TString name );

  TMVA::MethodBDT* BookMethod( Factory* factory , DataLoader* dataloader );
  

  double ROCIntegral;
  double Signal_TrainTest_Chi2 , Signal_TrainTest_Kolmo , Bkg_TrainTest_Kolmo , Bkg_TrainTest_Chi2 ; 
  void CalcEfficiency(TMVA::IMethod * method , TDirectory* dir , TString dsname );
  
  double ReadROC(TMVA::IMethod * method  );
  double ReadOverTrainingParam(TDirectory* dir , TString methodTitle , TString dsname , bool signal , bool kolmo );

  void PrintAll(std::ostream& out  , bool header = false){
    if( header )
      out << "name" << ","
	   << "ntrees" << ","
	   << "minnodesize" << ","
	   << "maxdepth" << ","
	   << "adaboostbeta" << ","
	   << "ncuts" << ","
	   << "ROCIntegral" << ","
	   << "Signal_TrainTest_Kolmo" << ","
	   << "Signal_TrainTest_Chi2" << ","
	   << "Bkg_TrainTest_Kolmo" << ","
	   << "Bkg_TrainTest_Chi2" << endl;
    else
      out << name << ","
	   << ntrees << ","
	   << minnodesize << ","
	   << maxdepth << ","
	   << adaboostbeta << ","
	   << ncuts << ","
	   << ROCIntegral << ","
	   << Signal_TrainTest_Kolmo << ","
	   << Signal_TrainTest_Chi2 << ","
	   << Bkg_TrainTest_Kolmo << ","
	   << Bkg_TrainTest_Chi2 << endl;
  }
};

class BDTOptimizer : public std::vector<BDTOptions> {
public:
  TString dsName;
  TString MVATitle;
  
  BDTOptimizer( TString dsname , TString MVATitle_ , string options_ );
  void FillRange( vector<string> range_info , vector<int>& vals );
  void FillRange( vector<string> range_info , vector<double>& vals );

  void EvaluateAll( Factory* factory , TDirectory* dir ){
    if( this->size() == 0 )
      return;
    
    ofstream myfile;
    myfile.open ( MVATitle + ".csv" );
    this->at(0).PrintAll(myfile , true);
    for(auto option : *this){
      option.CalcEfficiency( factory->GetMethod( dsName, option.name ) , dir , dsName ) ;
      option.PrintAll(myfile);
    }
  }
};

const string DEFAULT_INFNAME  = "/home/hbakhshi/Desktop/tHq/nTuples/FoxWolfram2";

int TMVAClassification( TString myMethodList , TString extention, BDTOptimizer* tth_bdt_options , BDTOptimizer* ttgj_bdt_options, BDTOptimizer* dig_bdt_options  )
{
  TMVA::MethodBDT* method_ttH ,* method_DiG ,* method_ttGj;
  
   TMVA::Tools::Instance();
   std::map<std::string,int> Use;
   Use["ttH"] = 1;
   Use["DiG"] = 1;
   Use["ttGj"] = 1;
   
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++)
	it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
	std::string regMethod(mlist[i]);

	if (Use.find(regMethod) == Use.end()) {
	  std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
	  for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
	  std::cout << std::endl;
	  return 1;
	}
	Use[regMethod] = 1;

	// if(regMethod == "ttH")
	//   extention += "_" + tth_bdt_option.to_string("tth");
	// else if(regMethod == "DiG")
	//   extention += "_" + dig_bdt_option.to_string("dig");
	// else if(regMethod == "ttGj")
	//   extention += "_" + ttgj_bdt_option.to_string("ttgj");
      }
   }

   TString infname     = DEFAULT_INFNAME;

   TString outfileName( "TMVA_"+ extention  +".root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   
   TMVA::Factory *factory = new TMVA::Factory( extention , outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=D:AnalysisType=Classification" );
   
   TFile *inputS = TFile::Open( infname + "/Signal.root" );
   TTree *signalTree     = (TTree*)inputS->Get("tHq/Trees/Events");
   Double_t signalWeight     = 1.0;
   TString default_w_str = "Weight.W0 * G1.w * G2.w";
   
   TMVA::DataLoader *dataloader_tth = NULL;
   TMVA::DataLoader *dataloader_dig = NULL;
   TMVA::DataLoader *dataloader_ttgg = NULL;
   
   if( Use["ttH"] ){

     std::cout << "Loading ttH trees" << endl;
     TFile *inputB = TFile::Open( infname + "/ttH.root" );
     TTree *background_ttH  = (TTree*)inputB->Get("tHq/Trees/Events");

     dataloader_tth =new TMVA::DataLoader(tth_bdt_options->dsName);
     // (please check "src/Config.h" to see all available global options)
     //
     //(TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
     //(TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

     dataloader_tth->AddVariable( "nJets", "nJets", "", 'I' ) ;
     dataloader_tth->AddVariable( "Max$( abs(jetsEta) )","jprimeeta" , "" , 'F' );
     dataloader_tth->AddVariable( "met.pt", "met", "", 'F' );
     dataloader_tth->AddVariable( "lepton.charge","LepCharge" , "" , 'I' );
     dataloader_tth->AddVariable( "eventshapes.aplanarity" , "aplanarity" , "" , 'F' );
     dataloader_tth->AddVariable( "foxwolf1.ONE" , "fwf1ONE" , "" , 'F' );
     
     // dataloader->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
     // dataloader->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );
     
     dataloader_tth->AddSignalTree( signalTree,     signalWeight );
     dataloader_tth->AddBackgroundTree( background_ttH, 1 );

     dataloader_tth->SetBackgroundWeightExpression( default_w_str  );
     dataloader_tth->SetSignalWeightExpression( default_w_str );
     TString cut_tth = "(DiG.mass > 100) && (Sum$(jetsPt>30) > 1 ) && (nMbJets==1) && (jetsPt[0] > 30)  && (LeptonType == 1 || LeptonType == 2) && (lepton.pt > 20)";
     TCut mycuts_tth = TCut(cut_tth);
     TCut mycutb_tth = TCut(cut_tth);
     dataloader_tth->PrepareTrainingAndTestTree( mycuts_tth, mycutb_tth,
						 "nTrain_Signal=3000:nTrain_Background=5000:SplitMode=Random:NormMode=NumEvents:!V");

     for(auto tth_bdt_option : *tth_bdt_options){
       method_ttH = tth_bdt_option.BookMethod( factory , dataloader_tth );
       cout << method_ttH->GetName() << endl;
     }
   }
   
   if( Use["DiG"] ){
     std::cout << "Loading DiG trees" << endl;
     TChain *backgroundDiG = new TChain("tHq/Trees/Events");
     backgroundDiG->Add(infname + "/DiPhoton_Jets.root" );


     dataloader_dig=new TMVA::DataLoader(dig_bdt_options->dsName);

     dataloader_dig->AddVariable( "nJets", "nJets", "", 'I' ) ;
     dataloader_dig->AddVariable( "met.pt", "met", "", 'F' );
     dataloader_dig->AddVariable( "DiG.pt" , "digpt" , "" , 'F' );
     dataloader_dig->AddVariable( "abs(DiG.eta)" , "digeta" , "" , 'F' );
     dataloader_dig->AddVariable( "DiG.mva" , "digmva" , "" , 'F' );
     dataloader_dig->AddVariable( "((G1.mva>G2.mva)*G2.mva +  (G1.mva<=G2.mva)*G1.mva)" , "minGmva" , "" , 'F' );
     dataloader_dig->AddSignalTree( signalTree,     signalWeight );
     dataloader_dig->AddBackgroundTree( backgroundDiG, 1 );

     dataloader_dig->SetBackgroundWeightExpression( default_w_str  );
     dataloader_dig->SetSignalWeightExpression( default_w_str );
     TString cut_dig = "(DiG.mass > 100)"; // && (met > 30)";
     TCut mycuts_dig = TCut(cut_dig);
     TCut mycutb_dig = TCut(cut_dig);
     dataloader_dig->PrepareTrainingAndTestTree( mycuts_dig, mycutb_dig,
						 "nTrain_Signal=10000:nTrain_Background=10000:SplitMode=Random:NormMode=NumEvents:!V");
     
     for(auto dig_bdt_option : *dig_bdt_options){
       method_DiG = dig_bdt_option.BookMethod( factory , dataloader_dig );
       cout << method_DiG->GetName() << endl;
     }
   }

   if( Use["ttGj"] ){
     std::cout << "Loading ttGj trees" << endl;
     TChain *backgroundTTGG = new TChain("tHq/Trees/Events");
     backgroundTTGG->Add( infname + "/TTGJ.root" );
     backgroundTTGG->Add( infname + "/TGJ.root" );

     dataloader_ttgg=new TMVA::DataLoader(ttgj_bdt_options->dsName);

     dataloader_ttgg->AddVariable( "nJets", "nJets", "", 'I' ) ;
     dataloader_ttgg->AddVariable( "met.pt", "met", "", 'F' );
     dataloader_ttgg->AddVariable( "DiG.pt" , "digpt" , "" , 'F' );
     dataloader_ttgg->AddVariable( "abs(DiG.eta)" , "digeta" , "" , 'F' );
     dataloader_ttgg->AddVariable( "DiG.mva" , "digmva" , "" , 'F' );
     dataloader_ttgg->AddVariable( "((G1.mva>G2.mva)*G2.mva +  (G1.mva<=G2.mva)*G1.mva)" , "minGmva" , "" , 'F' );
     
     dataloader_ttgg->AddSignalTree( signalTree,     signalWeight );
     dataloader_ttgg->AddBackgroundTree( backgroundTTGG, 1 );

     dataloader_ttgg->SetBackgroundWeightExpression( default_w_str  );
     dataloader_ttgg->SetSignalWeightExpression( default_w_str );
     TString cut_ttgg = "(DiG.mass > 100) && (nMbJets==1)" ;//&& (LeptonType == 1 || LeptonType == 2) && (lepton.pt > 20)"; && (jetsPt[0] > 30) && (met > 30) 
     TCut mycuts_ttgg = TCut(cut_ttgg);
     TCut mycutb_ttgg = TCut(cut_ttgg);
     dataloader_ttgg->PrepareTrainingAndTestTree( mycuts_ttgg, mycutb_ttgg,
						  "nTrain_Signal=10000:nTrain_Background=3000:SplitMode=Random:NormMode=NumEvents:!V");


     for(auto ttgj_bdt_option : *ttgj_bdt_options){
       method_ttGj = ttgj_bdt_option.BookMethod( factory , dataloader_ttgg );
       cout << method_ttGj->GetName() << endl;
     } 
   }
     

   
   // --------------------------------------------------------------------------------------------------
   //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
   // STILL EXPERIMENTAL and only implemented for BDT's !
   //
   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","FitGA");
   //
   // --------------------------------------------------------------------------------------------------

   // Now you can tell the factory to train, test, and evaluate the MVAs
   //
   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   tth_bdt_options->EvaluateAll( factory , outputFile );
   ttgj_bdt_options->EvaluateAll( factory , outputFile );
   dig_bdt_options->EvaluateAll( factory , outputFile );
   
   outputFile->Close();

   delete factory;
   return 0;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }

  return elems;
}

BDTOptions::BDTOptions(  TString name_ , int ntrees_  , double minnodesize_ , int maxdepth_  , double adaboostbeta_  , int ncuts_  ){
  name = name_;
  ntrees = ntrees_;
  ncuts = ncuts_;
  adaboostbeta = adaboostbeta_;
  maxdepth = maxdepth_;
  minnodesize = minnodesize_;
}

void BDTOptions::Set( string options_ ){
  vector<string> options =  split( options_ , ':' );
  for(auto opt_ : options){
    vector<string> opt =  split( opt_ , '=' ) ;
    if(opt.size() != 2 ){
      cout << opt_ << " is a wrong word" << endl;
      continue;
    }
    if( opt[0] == "nt" )
      ntrees = stoi( opt[1] );
    else if( opt[0] == "mns" )
      minnodesize = stod( opt[1] );
    else if( opt[0] == "md" )
      maxdepth = stoi( opt[1] );
    else if( opt[0] == "abb" )
      adaboostbeta = stod( opt[1] );
    else if( opt[0] == "nc" )
      ncuts = stoi( opt[1] );
    else
      cout << opt[0] << " is a wrong word " << endl;
  }
}
  
TString BDTOptions::to_string( TString name ){
  return TString( name + "_nt" + std::to_string(ntrees)
		  + "_mns" + TString::Format("%.1f",minnodesize)
		  + "_md" + std::to_string(maxdepth)
		  + "_abb" + TString::Format("%.1f" , adaboostbeta)
		  + "_nc" + std::to_string(ncuts) );
}


BDTOptimizer::BDTOptimizer( TString dsname , TString MVATitle_ , string options_ ){
  dsName = dsname ;
  MVATitle = MVATitle_;


  vector<int> ntrees ;
  vector<double> minnodesize ;
  vector<int> maxdepth ;
  vector<double> adaboostbeta ;
  vector<int> ncuts ;

  cout << options_ <<endl;
  
  vector<string> options =  split( options_ , ':' );
  for(auto opt_ : options){
    cout << opt_ <<endl;
    vector<string> opt =  split( opt_ , '=' ) ;
    if(opt.size() != 2 ){
      cout << opt_ << " is a wrong word" << endl;
      continue;
    }
    cout << opt[0] << opt[1] << endl;
    vector<string> range_info = {opt[1]};
    if(opt[1].find( ',' ) != string::npos)
      range_info = split( opt[1] , ',' );
      
    if( opt[0] == "nt" )
      FillRange( range_info , ntrees);
    else if( opt[0] == "mns" )
      FillRange( range_info , minnodesize);
    else if( opt[0] == "md" )
      FillRange( range_info , maxdepth);
    else if( opt[0] == "abb" )
      FillRange( range_info , adaboostbeta);
    else if( opt[0] == "nc" )
      FillRange( range_info , ncuts);
    else
      cout << opt[0] << " is a wrong word " << endl;

    cout << ntrees.size() << " " << minnodesize.size() << " " << maxdepth.size() << " " << adaboostbeta.size() << " " << ncuts.size() << endl;
  }

  for(auto nt : ntrees)
    for(auto mns : minnodesize)
      for(auto md : maxdepth)
	for(auto adabb : adaboostbeta)
	  for(auto nc : ncuts)
	    this->emplace_back( MVATitle + to_string( this->size() ) ,  nt , mns , md , adabb , nc );
	    
  for(auto i : *this)
    i.PrintAll(cout);
	    
}

void BDTOptimizer::FillRange( vector<string> range_info , vector<int>& vals ){
  if(range_info.size() == 1)
    vals.push_back( stoi(range_info[0]) );
  else if(range_info.size() != 3){
    vals.push_back(-1);
    cout << "wrong range info : " << range_info[0] << endl;
  }
  else{
    int from = stoi(range_info[0]);
    int to = stoi(range_info[1]);
    int step = stoi(range_info[2]);

    while( from <= to ){
      vals.push_back( from );
      from += step;
    }
  }
}
void BDTOptimizer::FillRange( vector<string> range_info , vector<double>& vals ){
  if(range_info.size() == 1)
    vals.push_back( stod(range_info[0]) );
  else if(range_info.size() != 3){
    vals.push_back(-1.);
    cout << "wrong range info : " << range_info[0] << endl;
  }
  else{
    double from = stod(range_info[0]);
    double to   = stod(range_info[1]);
    double step = stod(range_info[2]);

    while( from <= to ){
      vals.push_back( from );
      from += step;
    }
  }
}

void BDTOptions::CalcEfficiency(TMVA::IMethod * method , TDirectory* dir , TString dsname ){
  ROCIntegral = ReadROC( method );

  Signal_TrainTest_Kolmo = ReadOverTrainingParam(dir , method->GetName() , dsname , true , true );
  Signal_TrainTest_Chi2 = ReadOverTrainingParam(dir , method->GetName() , dsname , true , false );

  Bkg_TrainTest_Kolmo = ReadOverTrainingParam(dir , method->GetName() , dsname , false , true );
  Bkg_TrainTest_Chi2 = ReadOverTrainingParam(dir , method->GetName() , dsname , false , false );
}

double BDTOptions::ReadROC(TMVA::IMethod * method  ){
  TString title = method->GetName() ;
  MethodBase* theMethod = dynamic_cast<MethodBase*>(method);
  if(theMethod==0) return -1.0;
  TMVA::Results *results= theMethod->Data()->GetResults(name,Types::kTesting,Types::kClassification);
  std::vector<Float_t> *mvaRes = dynamic_cast<ResultsClassification *>(results)->GetValueVector();
  std::vector<Bool_t>  *mvaResType = dynamic_cast<ResultsClassification *>(results)->GetValueVectorTypes();
  Double_t fROCalcValue = 0;
  TMVA::ROCCurve *fROCCurve = nullptr;
  if (mvaResType->size() != 0) { 
    fROCCurve = new TMVA::ROCCurve(*mvaRes, *mvaResType);
    fROCalcValue = fROCCurve->GetROCIntegral();
  }
  return fROCalcValue ;           
}

double BDTOptions::ReadOverTrainingParam(TDirectory* dir , TString methodTitle , TString dsname , bool signal , bool kolmo ){
  TString hname = dsname + "/Method_" + methodTitle + "/" + methodTitle + "/MVA_" + methodTitle;
  TString appendix = signal ? "_S" : "_B" ;
  TH1* Test  = (TH1*) ( dir->Get( hname + appendix ) );
  TH1* Train = (TH1*) ( dir->Get( hname + "_Train" + appendix ) );
    
  TMVAGlob::NormalizeHists( Test , Train );

  for(int i = 0 ; i < Test->GetNbinsX() + 1 ; i ++){
    if( Test->GetBinContent(i) < 0 )
      Test->SetBinContent( i , 0 );
    if( Train->GetBinContent(i) < 0 )
      Train->SetBinContent( i , 0 );
  }
  
  Double_t ret = 0;
  if(kolmo )
    ret = Test->KolmogorovTest( Train, "X" ) ;
  else
    ret = Test->Chi2Test( Train , "WW CHI2/NDF" ); 
  return ret;
}

TMVA::MethodBDT* BDTOptions::BookMethod( Factory* factory , DataLoader* dataloader ){
  TMVA::MethodBase* method = factory->BookMethod( dataloader, TMVA::Types::kBDT, name,
						  "NegWeightTreatment=InverseBoostNegWeights:!H:!V:BoostType=AdaBoost:SeparationType=GiniIndex:CreateMVAPdfs:nCuts="+to_string(ncuts) );
  MethodBDT* ret = dynamic_cast<TMVA::MethodBDT*>( method );
  
  ret->SetNTrees(ntrees);
  ret->SetMinNodeSize(minnodesize);
  ret->SetMaxDepth(maxdepth);
  ret->SetAdaBoostBeta(adaboostbeta);
  
  return ret;
}



int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
  BDTOptimizer* tth_bdt_options,* ttgj_bdt_options , * dig_bdt_options;
  tth_bdt_options = ttgj_bdt_options = dig_bdt_options = NULL;
  TString methodList,ext;
  for (int i=1; i<argc; i++) {
    TString regMethod(argv[i]);
    if( regMethod=="--ext"){
	i++;
	ext = TString(argv[i]);
	continue;
    }
    if( regMethod=="--tth"){
      i++;
      tth_bdt_options = new BDTOptimizer("ttH" , "BDT_TTH" , string(argv[i] ) );
      for( auto i : *tth_bdt_options)
	i.PrintAll(cout);
      methodList += "ttH";
      continue;
    }
    if( regMethod=="--dig"){
      i++;
      dig_bdt_options = new BDTOptimizer("DiG" , "BDT_DiG" , string(argv[i] ) );
      methodList += "DiG";
      continue;
    }
    if( regMethod=="--ttgj"){
      i++;
      ttgj_bdt_options = new BDTOptimizer("ttGj" , "BDT_ttGj" , string(argv[i] ) );
      methodList += "ttGj";
      continue;
    }
  }
  
  return TMVAClassification(methodList , ext , tth_bdt_options , ttgj_bdt_options , dig_bdt_options);
}
