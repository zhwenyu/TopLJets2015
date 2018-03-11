#include "TopLJets2015/TopAnalysis/interface/KinematicsTools.h"
#include <iostream>


//
JetPullInfo_t getPullVector( MiniEvent_t &ev, int ijet)
{
  JetPullInfo_t result;
  result.n=0; result.nch=0;
  result.pull=TVector2(0,0);
  result.chPull=TVector2(0,0);
  
  //re-reconstruct the jet direction with the charged tracks
  TLorentzVector jet(0,0,0,0);
  jet.SetPtEtaPhiM(ev.j_pt[ijet], ev.j_eta[ijet], ev.j_phi[ijet], ev.j_mass[ijet]);
  TLorentzVector chargedJet(0,0,0,0);
  TLorentzVector constituent(0,0,0,0);
  std::vector<std::pair<TLorentzVector,bool> > allConstituents;
  unsigned int nCharged = 0;
  for(Int_t idx = 0; idx<ev.npf; idx++)
    {
      if(ev.pf_j[idx]!=ijet) continue;
      constituent.SetPtEtaPhiM( ev.pf_pt[idx], ev.pf_eta[idx], ev.pf_phi[idx], ev.pf_m[idx]);
      bool isCharged(abs(ev.pf_id[idx])==11 ||
		     abs(ev.pf_id[idx])==13 ||
		     abs(ev.pf_id[idx])==211 );
      allConstituents.push_back(std::make_pair(constituent,isCharged) );
      if(isCharged)
	{
	  chargedJet += constituent;
	  ++nCharged;      
	}
    }
  std::cout << std::endl;
  result.n=(Int_t) allConstituents.size();
  result.nch=nCharged;

  //stop here if <2 charged
  if( nCharged < 2 ) return result;

  //compute the pull
  double jetPt        = jet.Pt(),        jetPhi=jet.Phi(),                 jetRapidity=jet.Rapidity();
  double jetPtCharged = chargedJet.Pt(), jetPhiCharged = chargedJet.Phi(), jetRapidityCharged = chargedJet.Rapidity();
  TVector2 r(0,0);
  TVector2 pullAll(0,0);
  TVector2 pullCharged(0,0);
  for(size_t idx = 0; idx<allConstituents.size(); ++idx)
    {
      TLorentzVector &cp4=allConstituents[idx].first;
      bool &isCharged=allConstituents[idx].second;
      double constituentPt       = cp4.Pt();
      double constituentPhi      = cp4.Phi();
      double constituentRapidity = cp4.Rapidity();
      r.Set( constituentRapidity - jetRapidity, TVector2::Phi_mpi_pi( constituentPhi - jetPhi ) );
      pullAll += ( constituentPt / jetPt ) * r.Mod() * r;
      //calculate TVector using only charged tracks
      if( isCharged )
	r.Set( constituentRapidity - jetRapidityCharged, TVector2::Phi_mpi_pi( constituentPhi - jetPhiCharged ) );
      pullCharged += ( constituentPt / jetPtCharged ) * r.Mod() * r;
    }
  
  result.pull=pullAll;
  result.chPull=pullCharged;
  return result;
}


//
Float_t computeMT(TLorentzVector &a, TLorentzVector &b)
{
  return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
}

