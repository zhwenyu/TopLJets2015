import pickle
import ROOT
import random

class EventMixingTool:

    def __init__(self, mixFile):

        """ Reads the event mixing data from a pickle file and builds a list of crossing angle probabilities """

        self.mixedRP=None
        self.xangleRelFracs={}

        try:
            print '[EventMixingTool] analysing mixfile',mixFile
            with open(mixFile,'r') as cachefile:
                self.mixedRP=pickle.load(cachefile)        

            #build the list of probabilities for the crossing angles in each era
            for key in mixedRP:
                mix_era,mix_xangle,mix_evcat=key
                if not mix_xangle in VALIDLHCXANGLES: continue
                n=len(mixedRP[key])
                xangleKey=(mix_era,mix_evcat)
                if not xangleKey in xangleRelFracs:
                    xangleRelFracs[xangleKey]=ROOT.TH1F('xanglefrac_%s_%d'%xangleKey,'',len(VALIDLHCXANGLES),0,len(VALIDLHCXANGLES))
                xbin=(mix_xangle-120)/10+1
                xangleRelFracs[xangleKey].SetBinContent(xbin,n)

            for xangleKey in xangleRelFracs:
                xangleRelFracs[xangleKey].Scale(1./xangleRelFracs[xangleKey].Integral())
            
            print '[EventMixingTool] events will be mixed from',mixFile
        except Exception as e:
            if mixFile : print e
            pass


    def isIdle(self):
        
        """flag if it didn't initialize correctly"""

        return False if self.mixedRP else True


    def getNew(self,evEra,beamXangle,isData,validAngles,mixEvCategs):

        """get new list of mixed protons from different event categories
        also returns a list of variables which be used for pileup discrimination"""

        mixed_far_rptks={}
        mixed_near_rptks={}
        mixed_pudiscr={}

        try:
            for mixEvCat in mixEvCategs:
                mixed_far_rptks[mixEvCat]=None
                mixed_near_rptks[mixEvCat]=None
                mixed_pudiscr[mixEvCat]=None

                if isData and not beamXangle in validAngles : continue

                mixedEvKey                 = (evEra,beamXangle,mixEvCat)
                mixedEv                    = random.choice( self.mixedRP[mixedEvKey] )
                mixed_pudiscr[mixEvCat]    = mixedEv.puDiscr
                mixed_far_rptks[mixEvCat]  = mixedEv.far_rptks
                mixed_near_rptks[mixEvCat] = mixedEv.near_rptks
        except Exception as e:
            print e  
            print evEra,beamXangle,mixEvCat,'->',mixedEvKey
            pass


    def mergeWithMixedEvent(self,far_rptks,mixed_far_rptks,near_rptks,mixed_near_rptks):
            
        """merges tracks from two different events (useful for pure signal embedding) """
                    
        for mixEvCat in mixed_far_rptks:
            tksPos=mixed_far_rptks[mixEvCat][0]+far_rptks[0]
            shuffle(tksPos)
            tksNeg=mixed_far_rptks[mixEvCat][1]+far_rptks[1]
            shuffle(tksNeg)
            mixed_far_rptks[mixEvCat]=(tksPos,tksNeg)

            tksPos=mixed_near_rptks[mixEvCat][0]+near_rptks[0]
            shuffle(tksPos)
            tksNeg=mixed_near_rptks[mixEvCat][1]+near_rptks[1]
            shuffle(tksNeg)
            mixed_near_rptks[mixEvCat]=(tksPos,tksNeg)
            
        return mixed_far_rptks,mixed_near_rptks

