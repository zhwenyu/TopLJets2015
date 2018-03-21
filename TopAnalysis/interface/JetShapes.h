#ifndef _jetshapes_h_
#define _jetshapes_h_

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

double calcGA(double beta, double kappa, pat::Jet *jet, bool includeNeutrals, double ptcut);
double getMult(pat::Jet *jet, bool includeNeutrals, double ptcut);
double getPtD(pat::Jet *jet, bool includeNeutrals, double ptcut);
double getPtDs(pat::Jet *jet, bool includeNeutrals, double ptcut);
double getWidth(pat::Jet *jet, bool includeNeutrals, double ptcut);
double getEcc(pat::Jet *jet, bool includeNeutrals, double ptcut);
double getTau(int N, int M, pat::Jet *jet, bool includeNeutrals, double ptcut);
double getC(int N, double beta, pat::Jet *jet, bool includeNeutrals, double ptcut);
std::vector<double> getZg(pat::Jet *jet, bool includeNeutrals, double ptcut);
double mapAngleMPiToPi(double phi);
double getNSD(double zcut, double beta, pat::Jet *jet, bool includeNeutrals, double ptcut);
std::map<TString,double> getECF(pat::Jet *jet, bool includeNeutrals, double ptcut);
double getPFFraction(std::vector<int> pids, pat::Jet *jet);
double mapAngleMPiToPi(double phi);

#endif
