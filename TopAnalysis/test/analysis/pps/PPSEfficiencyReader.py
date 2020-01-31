import ROOT
import sys

class PPSEfficiencyReader:
    
    """ 
    takes care of reading the efficency measurements to memory and retrieving the final efficiency correction 
    see details in https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsStripsEfficiencies
    """

    def __init__(self, fList, year=2017):

        self.allEffs={}

        for fIn in fList.split(','):
            baseDir='Strips/%d'%year if 'MultiTrack' in fIn else 'Pixel/%d'%year
            fIn=ROOT.TFile.Open(ROOT.gSystem.ExpandPathName(fIn))
            for k in fIn.Get(baseDir).GetListOfKeys():
                for kk in fIn.Get(baseDir+'/'+k.GetName()).GetListOfKeys():
                    hname=kk.GetName()
                    if '2D'in hname : continue
                    self.allEffs[hname]=kk.ReadObj()
                    self.allEffs[hname].SetDirectory(0)
            fIn.Close()

        print '[PPSEfficiencyReader] retrieved %d histograms'%len(self.allEffs)


    def getPPSEfficiency(self,era,xangle,xi,rp,isMulti=True, applyMultiTrack=False):

        sector=45 if rp<100 else 56
        eff,effUnc=1.0,0.0
        if isMulti:

            if applyMultiTrack:
                multiTrack=self.allEffs['h%dmultitrackeff_%s_avg_RP%d'%(sector,era,rp)]
                ieff = multiTrack.GetBinContent(1)
                eff *= ieff
        
            raddam=self.allEffs['h%d_%s_%d_1D'%(sector,era,xangle)]
            raddamUnc=self.allEffs['h%derrors_%s_%d_1D'%(sector,era,xangle)]
            ibin=raddam.FindBin(xi)
            ieff    = raddam.GetBinContent(ibin)
            if ieff>0:
                eff    *= ieff
                effUnc += (raddamUnc.GetBinError(ibin)/ieff)**2            
                
        else:

            # FIXME
            # pixels are not  fully available yet so assume 100% eff
            eff,effUnc=1.0,0.0

            
        effUnc=eff*ROOT.TMath.Sqrt(effUnc)

        return eff,effUnc


    def getCombinedEfficiency(self,ppsPosEff,ppsPosEffUnc,ppsNegEff,ppsNegEffUnc,i_proton_cat):
        ppsEff,ppsEffUnc=0.,0.
        if i_proton_cat==1:
            ppsEff=ppsPosEff*ppsNegEff
            ppsEffUnc=ROOT.TMath.Sqrt( (ppsPosEff*ppsNegEffUnc)**2+(ppsPosEffUnc*ppsNegEff)**2)
        elif i_proton_cat==2:
            ppsEff=ppsPosEff
            ppsEffUnc=ppsPosEffUnc
        elif i_proton_cat==3:
            ppsEff=ppsNegEff
            ppsEffUnc=ppsNegEffUnc
        return ppsEff,ppsEffUnc


def main():

    ppEffReader=PPSEfficiencyReader(fList='test/analysis/pps/PreliminaryEfficiencies_October92019_1D2DMultiTrack.root')

    #test for different conditions
    for era in ['2017B','2017C','2017D','2017E','2017F']:
        for xangle in [120,130,140,150]:
            for rp in [3,103]:
                eff=ppEffReader.getPPSEfficiency(era,xangle,0.035,rp)
                print '%6s %d %3d %3.3f +/- %3.3f'%(era,xangle,rp,eff[0],eff[1])


if __name__ == "__main__":
    sys.exit(main())
