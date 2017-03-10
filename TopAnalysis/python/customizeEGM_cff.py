import FWCore.ParameterSet.Config as cms
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

# EGM corrections
# scale regression https://twiki.cern.ch/twiki/bin/view/CMS/EGMRegression
# smearer https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer
def customizeEGM(process,runOnData):

    #smearing
    process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
    process.calibratedPatElectrons.correctionFile = 'EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_23Jan_ele'
    process.calibratedPatElectrons.isMC=False if runOnData else True
    print 'Using smeared electrons with corrections from',process.calibratedPatElectrons.correctionFile

    # Set up electron ID (VID framework)
    dataFormat = DataFormat.MiniAOD
    switchOnVIDElectronIdProducer(process, dataFormat)
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff'] 
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

    #if regression is applied vid must be fixed https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1837.html
    process = regressionWeights(process)
    process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')
    process.load('Configuration.StandardSequences.Services_cff')
    print 'EGM regression will be applied prior to smearing'
    process.electronRegressionValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
    process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
    process.egmSeq=cms.Sequence(process.regressionApplication*process.calibratedPatElectrons*process.egmGsfElectronIDSequence)
    #process.egmSeq=cms.Sequence(process.regressionApplication*process.egmGsfElectronIDSequence)
