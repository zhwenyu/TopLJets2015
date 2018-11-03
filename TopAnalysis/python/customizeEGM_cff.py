import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

# EGM corrections :  https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2

def customizeEGM(process,era):

    #switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
    #my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff', 
    #                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff', ]
    #for idmod in my_id_modules:
    #    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


    if not 'era2017' in era: return None
    setupEgammaPostRecoSeq(process,
                           runVID=True,
                           era='2017-Nov17ReReco')

    process.egammaPostReco=cms.Path(process.egammaPostRecoSeq)
    return process.egammaPostReco
