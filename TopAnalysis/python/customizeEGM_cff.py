import FWCore.ParameterSet.Config as cms
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

# EGM corrections
# scale regression https://twiki.cern.ch/twiki/bin/view/CMS/EGMRegression
# smearer https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer
def customizeEGM(process,runOnData):

    process = regressionWeights(process)
    process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')
    process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
    process.calibratedPatElectrons.isMC = cms.bool(False) if runOnData else cms.bool(True)
    print 'Using smeared electrons with corrections from',process.calibratedPatElectrons.correctionFile

    # Set up electron ID (VID framework)
    dataFormat = DataFormat.MiniAOD
    switchOnVIDElectronIdProducer(process, dataFormat)
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff']
    #'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff'] 
    for idmod in my_id_modules:
        print idmod
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    process.selectedElectrons = cms.EDFilter("PATElectronSelector",
                                             src = cms.InputTag("calibratedPatElectrons"),
                                             cut = cms.string("pt>5")
                                             )
    process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedElectrons')
    #process.electronIDValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')
    process.electronRegressionValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')
    process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')
    print 'VID setup to be applied on ',process.electronMVAValueMapProducer.srcMiniAOD

    process.EGMRegression = cms.Path(process.regressionApplication)
    process.EGMSmearerElectrons = cms.Path(process.calibratedPatElectrons)

    #process.egmSeq=cms.Sequence(process.regressionApplication*
    #                            process.calibratedPatElectrons*
    #                            process.selectedElectrons * 
    #                            process.egmGsfElectronIDSequence)
