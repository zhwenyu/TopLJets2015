import ROOT
import sys

class ParametrizedLeptonScaleUnc:
    """A wrapper to breakup lepton energy scale uncertainties for TOP-17-010"""

    def __init__(self):
        self.kamuca=ROOT.KalmanMuonCalibrator("MC_80X_13TeV")
        
        self.euncs=[]
        eurl='${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/era2016/ParametrizedLeptonScaleUnc.root'
        fIn=ROOT.TFile.Open(eurl)
        for i in xrange(0,3):
            self.euncs.append(fIn.Get('scaleUnc_11_%d'%i))
            self.euncs[-1].SetDirectory(0)
        fIn.Close()

    def getUncertainty(self,pdgId,basePt,baseEta,basePhi,charge):
        if abs(pdgId)==11:
            return self.getElectronUncertainty(basePt,baseEta)
        else:
            return self.getMuonUncertainty(basePt,baseEta,basePhi,charge)

    def getElectronUncertainty(self,basePt,baseEta):
        ptToEval=max(self.euncs[0].GetXaxis().GetXmin(),min(self.euncs[0].GetXaxis().GetXmax()-0.01,basePt))
        etaToEval=max(self.euncs[0].GetYaxis().GetXmin(),min(self.euncs[0].GetYaxis().GetXmax()-0.01,abs(baseEta)))
        relUncs=[ h.GetBinContent( h.GetXaxis().FindBin(ptToEval), h.GetYaxis().FindBin(etaToEval) )
                  for h in self.euncs ]
        return relUncs

    def getMuonUncertainty(self,basePt,baseEta,basePhi,charge):
        pt=self.kamuca.getCorrectedPt(basePt,baseEta,basePhi,charge)

        N=self.kamuca.getN()
        syst=0
        for i in range(0,N):
            self.kamuca.vary(i,+1)
            syst += (pt-self.kamuca.getCorrectedPt(basePt,baseEta,basePhi,charge))**2
        syst=ROOT.TMath.Sqrt(syst)/basePt

        self.kamuca.reset()
        self.kamuca.varyClosure(+1)
        syst_closure=abs(pt-self.kamuca.getCorrectedPt(basePt,baseEta,basePhi,charge))/basePt

        return [syst,syst_closure]

def main():
    
    plsu=ParametrizedLeptonScaleUnc()
    for args in [(11,40.,1.2,0,+1),
                 (13,40.,1.2,0,+1),
                 (13,40.,1.2,0,-1)]:
        print args
        print plsu.getUncertainty(*args)

if __name__ == "__main__":
    sys.exit(main())
