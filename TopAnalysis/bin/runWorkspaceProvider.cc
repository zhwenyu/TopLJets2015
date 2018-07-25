#include "TopLJets2015/TopAnalysis/interface/WorkspaceProvider.h"

using namespace std;

int main( int argc, char** argv )
{
  TString channel, histname, boson;
  YieldsErr YE;
  double sigEff =1;
  double bkgEff =1;
  int nBin = 4;
  for (int i=1; i<argc; i++) {
    TString input(argv[i]);
    if ( input=="--Chan"){
      i++;
      channel = TString(argv[i]);
      continue;
    } else if ( input=="--V"){
      i++;
      boson = TString(argv[i]);
      continue;
    } else if( input=="--Hist"){
      i++;
      histname = TString(argv[i]);
      continue;
    } else if( input=="--YieldErr"){
      i++;
      YE = splitter(string(argv[i]),':',',');     
      continue;
    } else if( input=="--sigEff"){
      i++;
      sigEff = (double)atof(TString(argv[i]));
      continue;
    } else if( input=="--bkgEff"){
      i++;
      bkgEff = (double)atof(TString(argv[i]));
      continue;
    } else if( input=="--nBin"){
      i++;
      nBin = (int)atof(TString(argv[i]));
      continue;
    }
  }
  
  WorkspaceProvider wsp(histname, channel, boson, nBin);
  wsp.ProvideWS();
  wsp.creatTFHists();
  wsp.makeCard(YE, sigEff, bkgEff);
  return 0;
}
