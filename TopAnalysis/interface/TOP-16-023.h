#ifndef _top16023_h_
#define _top16023_h_

#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"

void RunTop16023(TString filename,
		     TString outname,
		     Int_t channelSelection, 
		     Int_t chargeSelection, 
		     FlavourSplitting flavourSplitting,
		     TH1F *normH, 
		     Bool_t runSysts,
		     TString era);

#endif
