import ROOT
import sys
import numpy as np

def isPixelFiducial(era,sector,x,tx,y,ty):

    """
    check if the track is in the fiducial region
    cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/TaggedProtonsPixelEfficiencies
    """

    #check angle of the track to be below 20mrad
    if abs(tx)>0.02 : return False
    if abs(ty)>0.02 : return False         
        
    xy_fid=None
    if era in ['2017B','2017C','2017D']:
        if sector==45:
            xy_fid=[1.860,24.334,-11.098,4.298]
            if era=='2017B' :
                xy_fid[0]=1.995
                xy_fid[1]=24.479
        else:
            xy_fid=[2.422,24.620,-10.698,4.698]
    else:
        if sector==45:
            xy_fid=[1.995,24.479,-10.098,4.998]
        else:
            xy_fid=[2.422,24.620,-9.698,5.498]
    px_x0_rotated = x * np.cos((-8. / 180.) * np.pi) - y * np.sin((-8. / 180.) * np.pi)
    px_y0_rotated = x * np.sin((-8. / 180.) * np.pi) + y * np.cos((-8. / 180.) * np.pi)
    if px_x0_rotated<xy_fid[0] : return False
    if px_x0_rotated>xy_fid[1] : return False
    if px_y0_rotated<xy_fid[2] : return False
    if px_y0_rotated>xy_fid[3] : return False
    
    return True

def vetoPixels2017(run,era):

    #veto 'C2', 'D', 'F2', and 'F3' run ranges available in 
    #https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsPixelEfficiencies

    if 'D' in era: return True
    if run<0:
        #use random number assignment to throw away events in MC
        subEra=np.random.uniform()
        if 'C' in era and subEra>0.607 : return True
        if 'F' in era and subEra>0.128 : return True
    else:
        if 'C' in era and run>=300806 and run<=302029 : return True
        if 'F' in era and run>=305178:                  return True
    
    return False


def doFinalCheck2017(pos_protons,neg_protons,run,era):
    hasVeto=False
    if vetoPixels2017(run,era):
        if len(pos_protons[0])==0 and len(pos_protons[1])>0:
            pos_protons[1]=[]
            hasVeto=True
        if len(neg_protons[0])==0 and len(neg_protons[1])>0:
            neg_protons[1]=[]
            hasVeto=True
    return hasVeto,pos_protons,neg_protons


class PPSEfficiencyReader:
    
    """ 
    takes care of reading the efficency measurements to memory and retrieving the final efficiency correction 
    see details in https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsStripsEfficiencies
    """

    def __init__(self, fList, year=2017):

        #fix the seed for reproducible jobs
        np.random.seed(1)

        self.allEffs={}

        for fIn in fList.split(','):

            if 'MultiTrack' in fIn:
                baseDir='Strips/%d'%year
                fIn=ROOT.TFile.Open(ROOT.gSystem.ExpandPathName(fIn))
                for k in fIn.Get(baseDir).GetListOfKeys():
                    for kk in fIn.Get(baseDir+'/'+k.GetName()).GetListOfKeys():
                        hname=kk.GetName()
                        if '2D' in hname : continue
                        self.allEffs[hname]=kk.ReadObj()
                        self.allEffs[hname].SetDirectory(0)
                fIn.Close()
                
            elif 'pixelEfficiencies_multiRP' in fIn:

                fIn=ROOT.TFile.Open(ROOT.gSystem.ExpandPathName(fIn))
                for era in ['B','C1','C2','D','E','F1','F2','F3']:
                    baseDir='Pixel/{0}/{0}{1}'.format(year,era)
                    eraDir=fIn.Get(baseDir)
                    try:

                        for k in eraDir.GetListOfKeys():
                            hname=k.GetName()
                            if not '2D' in hname : continue
                            self.allEffs[hname+"_ip"]=k.ReadObj()
                            self.allEffs[hname+"_ip"].SetDirectory(0)
                            self.allEffs[hname+"_ip"].SetName(hname+"_ip")
                    except:
                        print era,'not available as',baseDir

                #compose lumi averaged for eras C and F
                for era,eras in [ ('C',[('C1',0.62),('C2',0.38)]),
                                  ('F',[('F1',0.13),('F2',0.59),('F3',0.28)]) ]:
                    
                    for rp in [45,56]:
                        hname='h{0}_220_2017{1}_all_2D_ip'
                        firstSubEra=eras[0][0]
                        inc_hname=hname.format(rp,era)
                        self.allEffs[inc_hname]=self.allEffs[hname.format(rp,firstSubEra)].Clone(inc_hname)
                        self.allEffs[inc_hname].Reset('ICE')
                        self.allEffs[inc_hname].SetDirectory(0)

                        for subEra,subEraWgt in eras:
                            self.allEffs[inc_hname].Add(self.allEffs[hname.format(rp,subEra)],subEraWgt)
                        
                fIn.Close()

            elif 'pixelEfficiencies_radiation' in fIn:

                #pixel efficiencies from radiation (not all eras available due to 3+3 readout)

                fIn=ROOT.TFile.Open(ROOT.gSystem.ExpandPathName(fIn))
                for era1, era2 in [('B','B'),('C','C1'),('E','E'),('F','F1')]:
                    for sector in [45,56]:
                        hname='h{0}_220_2017{1}_all_2D_pxrad'.format(sector,era)
                        self.allEffs[hname]=fIn.Get('Pixel/2017/2017{0}/h{1}_220_2017{0}_all_2D'.format(era2,sector))
                        self.allEffs[hname].SetDirectory(0)
                        self.allEffs[hname].SetName(hname+"_pxrad")
                fIn.Close()
        
        #pure 0 strip tracks eff from J. Kaspar
        e1f=7519./(7519.+1440.)
        self.pure0Probs={
            (45,120,'2017B'):0.8605,
            (45,120,'2017C'):0.8687,
            (45,120,'2017D'):0.8665,
            (45,120,'2017E'):e1f*1.0+(1-e1f)*0.6945,
            (45,120,'2017F'):0.6803,
            (45,130,'2017B'):0.7749,
            (45,130,'2017C'):0.7888,
            (45,130,'2017D'):0.7920,
            (45,130,'2017E'):e1f*1.0+(1-e1f)*0.4680,
            (45,130,'2017F'):0.4667,
            (45,140,'2017B'):0.7137,
            (45,140,'2017C'):0.7181,
            (45,140,'2017D'):0.7353,
            (45,140,'2017E'):e1f*1.0+(1-e1f)*0.3556,
            (45,140,'2017F'):0.3878,
            (45,150,'2017B'):0.6359,
            (45,150,'2017C'):0.6510,
            (45,150,'2017D'):0.6713,
            (45,150,'2017E'):e1f*1.0+(1-e1f)*0.3493,
            (45,150,'2017F'):0.3593,
            (56,120,'2017B'):0.8412,
            (56,120,'2017C'):0.8370,
            (56,120,'2017D'):0.8273,
            (56,120,'2017E'):e1f*0.6572+(1-e1f)*0.6307,
            (56,120,'2017F'):0.6053,
            (56,130,'2017B'):0.7409,
            (56,130,'2017C'):0.7400,
            (56,130,'2017D'):0.7375,
            (56,130,'2017E'):e1f*0.4822+(1-e1f)*0.3976,
            (56,130,'2017F'):0.3813,
            (56,140,'2017B'):0.6752,
            (56,140,'2017C'):0.6607,
            (56,140,'2017D'):0.6729,
            (56,140,'2017E'):e1f*0.3791+(1-e1f)*0.2982,
            (56,140,'2017F'):0.3100,
            (56,150,'2017B'):0.5948,
            (56,150,'2017C'):0.5896,
            (56,150,'2017D'):0.6010,
            (56,150,'2017E'):e1f*0.3467+(1-e1f)*0.2904,
            (56,150,'2017F'):0.2862,
            }

        print '[PPSEfficiencyReader] retrieved %d histograms'%len(self.allEffs)


    def getPPSEfficiency(self,era,xangle,xi,x,y,rp,isMulti=True, applyMultiTrack=False, applyInterPotAndPure0=True):

        sector=45 if rp<100 else 56
        eff,effUnc=1.0,0.0
        if isMulti:

            if applyMultiTrack:
                multiTrack=self.allEffs['h%dmultitrackeff_%s_avg_RP%d'%(sector,era,rp)]
                ieff = multiTrack.GetBinContent(1)
                eff *= ieff
                if ieff>0:
                    effUnc += (multiTrack.GetBinError(1)/ieff)**2

            raddam=self.allEffs['h%d_%s_%d_1D'%(sector,era,xangle)]
            raddamUnc=self.allEffs['h%derrors_%s_%d_1D'%(sector,era,xangle)]
            ibin=raddam.FindBin(xi)
            ieff = raddam.GetBinContent(ibin)
            eff *= ieff
            if ieff>0:
                effUnc += (raddamUnc.GetBinError(ibin)/ieff)**2            


            if applyInterPotAndPure0:

                pure0Eff = self.pure0Probs[(sector,xangle,era)]
                eff *= pure0Eff
                if pure0Eff>0:
                    effUnc *= 0.02**2 #ad-hoc

                if x>-90 and y>-90 : #-99 is the default for n/a
                    interPot=self.allEffs['h{0}_220_{1}_all_2D_ip'.format(sector,era)]
                    xbin=interPot.GetXaxis().FindBin(x)
                    ybin=interPot.GetYaxis().FindBin(y)
                    ieff = interPot.GetBinContent(xbin,ybin)
                    eff *=ieff
                    #the errors in the histograms seem flawed, assign 2% ad-hoc
                    effUnc += 0.02**2
                    #if ieff>0:
                    #    effUnc += (interPot.GetBinError(xbin,ybin)/ieff)**2

        else:

            #check if era is available (if not do nothing)
            hname='h{0}_220_{1}_all_2D_ip'.format(sector,era)       
            if hname in self.allEffs and  x>-90 and y>-90 :
                pxrad=self.allEffs[hname]
                xbin=pxrad.GetXaxis().FindBin(x)
                ybin=pxrad.GetYaxis().FindBin(y)
                ieff=pxrad.GetBinContent(xbin,ybin)
                eff *=ieff
                #the error in the histogram seems flawed, assign 2% ad-hoc
                effUnc += 0.02**2
                #if ieff>0:
                #    effUnc += (pxrad.GetBinError(xbin,ybin)/ieff)**2
            
        effUnc=eff*ROOT.TMath.Sqrt(effUnc)

        return eff,effUnc


    def getProjectedFinalState(self,
                               pos_protons,multiPosEff,multiPosEffUnc,pixelPosEff,pixelPosEffUnc,
                               neg_protons,multiNegEff,multiNegEffUnc,pixelNegEff,pixelNegEffUnc,
                               sighyp,run,era):

        """
        it receives a candidate list of protons for the positive and negative side
        and the corresponding efficiencies for multi- and pixel-only algorithms
        sighyp is a number between 0 and 16 where the bits represent
        0b - number of pixels in negative side
        1b - number of multi in negative side
        2b - number of pixels in positive side
        3b - number of multi in positive side
        it returns the final list of protons in each side corresponding to the sighyp
        and the final combine efficiency and its uncertainty
        """

        #check how many multi are required for the signal hypothesis
        nPixNegInSigHyp   = ((sighyp>>0) & 0x1)
        nMultiNegInSigHyp = ((sighyp>>1) & 0x1)
        nPixPosInSigHyp   = ((sighyp>>2) & 0x1)
        nMultiPosInSigHyp = ((sighyp>>3) & 0x1)

        #check how many multi and pixels are available
        nPixNeg   = min(1,len(neg_protons[1]))
        nMultiNeg = min(1,len(neg_protons[0]))
        nPixPos   = min(1,len(pos_protons[1]))
        nMultiPos = min(1,len(pos_protons[0]))

        #assign final proton multiplicities and compute efficiencies per arm
        pos_ppsWgt,pos_ppsWgtUnc,pos_protons = self.assignFinalSignalHypothesisToArm(nMultiPosInSigHyp, nMultiPos, multiPosEff, multiPosEffUnc,
                                                                                     nPixPosInSigHyp,   nPixPos,   pixelPosEff, pixelPosEffUnc,
                                                                                     pos_protons)

        neg_ppsWgt,neg_ppsWgtUnc,neg_protons = self.assignFinalSignalHypothesisToArm(nMultiNegInSigHyp, nMultiNeg, multiNegEff, multiNegEffUnc,
                                                                                     nPixNegInSigHyp,   nPixNeg,   pixelNegEff, pixelNegEffUnc,
                                                                                     neg_protons)
        
        #global PPS efficienciy
        ppsWgt    = pos_ppsWgt*neg_ppsWgt
        ppsWgtUnc = ROOT.TMath.Sqrt(pos_ppsWgtUnc**2+neg_ppsWgtUnc**2)

        #final check, if only pixels remain, ensure the run/era is not to be vetoed
        hasVeto,pos_protons,neg_protons = doFinalCheck2017(pos_protons,neg_protons,run,era)
        if hasVeto: 
            ppsWgt,ppsWgtUnc=0.0,0.0

        #return final result
        return pos_protons,neg_protons,ppsWgt,ppsWgtUnc

    def assignFinalSignalHypothesisToArm(self,
                                         nMultiInSigHyp, nMulti, multiEff, multiEffUnc,
                                         nPixInSigHyp,   nPix,   pixelEff, pixelEffUnc,
                                         protons):

        """ applies the final signal reconstruction assignments and probabilities for one arm """

        #kill protons in the list if needed
        if nMultiInSigHyp==0:
            protons[0]=[]
            protons[2]=[]
        if nPixInSigHyp==0:
            protons[1]=[]

        #compute arm weight (probability)
        ppsArmWgt,ppsArmWgtUnc = 0.0, 0.0

        #reconstructed |1 1>
        if nMulti==1 and nPix==1:
            
            #target |1 1>
            if nMultiInSigHyp==1 and nPixInSigHyp==1:
                ppsArmWgt    = multiEff*pixelEff
                ppsArmWgtUnc = (multiEffUnc/multiEff)**2 + (pixelEffUnc/pixelEff)**2

            #target |0 1>
            if nMultiInSigHyp==0 and nPixInSigHyp==1:
                ppsArmWgt    = (1-multiEff)*pixelEff
                ppsArmWgtUnc = (multiEffUnc/(1-multiEff))**2 + (pixelEffUnc/pixelEff)**2

            #target |0 0>
            if nMultiInSigHyp==0 and nPixInSigHyp==0:
                ppsArmWgt    = (1-multiEff)*(1-pixelEff)
                ppsArmWgtUnc = (multiEffUnc/(1-multiEff))**2 + (pixelEffUnc/(1-pixelEff))**2


        #reconstructed |0 1>
        if nMulti==0 and nPix==1:

            #target |0 1>
            if nMultiInSigHyp==0 and nPixInSigHyp==1:
                ppsArmWgt    = pixelEff
                ppsArmWgtUnc = (pixelEffUnc/pixelEff)**2

            #target |0 0>
            if nMultiInSigHyp==0 and nPixInSigHyp==0:
                ppsArmWgt    = (1-pixelEff)
                ppsArmWgtUnc = (pixelEffUnc/(1-pixelEff))**2

        #reconstructed |0 0>
        if nMulti==0 and nPix==0:

            #target |0 0>
            if nMultiInSigHyp==0 and nPixInSigHyp==0:
                ppsArmWgt    = 1.
                ppsArmWgtUnc = 0.

        #finalize error
        ppsArmWtUnc = ppsArmWgt*ROOT.TMath.Sqrt(ppsArmWgtUnc)

        return ppsArmWgt,ppsArmWgtUnc,protons


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
