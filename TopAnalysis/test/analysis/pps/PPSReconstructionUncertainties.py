import ROOT
import sys
import os

class PPSReconstructionUncertainty:

    """wraps up the procedure to evaluate the xi reconstruction uncertainties"""
    def __init__(self,era):

        self.uncs={}

        fIn=ROOT.TFile.Open('{}/src/TopLJets2015/TopAnalysis/test/analysis/pps/reco_charactersitics_version1.root'.format(os.environ['CMSSW_BASE']))
        
        for k in fIn.GetListOfKeys():

            period=k.GetName()
            if not str(era) in period: continue

            for kk in k.ReadObj().GetListOfKeys():
                rec=kk.GetName()
                if not rec in self.uncs: self.uncs[rec]={}
                if not period in self.uncs[rec] : self.uncs[rec][period]={}

                for gname in ['bias','resolution','systematics']:
                    self.uncs[rec][period][k]=fIn.Get('{}/{}/xi/g_{}_vs_xi'.format(period,rec,gname))

        fIn.Close()



def main():

    ppRecUnc=PPSReconstructionUncertainty(2017)


if __name__ == "__main__":
    sys.exit(main())
