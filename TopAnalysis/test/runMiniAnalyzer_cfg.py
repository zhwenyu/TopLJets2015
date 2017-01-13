import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnData', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run this on real data"
                 )
options.register('era', 'era2016',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "era to use (configurations defined in python/EraConfig.py)"
                 )
options.register('outFilename', 'MiniEvents.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('baseJetCollection','slimmedJetsPuppi',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Base jet collection"
                 )
options.register('inputDir', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input directory with files to process"
                 )
options.register('lumiJson', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'apply this lumi json file'
                 )
options.register('saveTree', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "save summary tree"
                 )
options.register('savePF', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'save PF candidates'
                 )
options.register('applyLeptonCorrections', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'Switch to apply lepton corrections'
                 )
options.parseArguments()

#get the configuration to apply
from TopLJets2015.TopAnalysis.EraConfig import getEraConfiguration
globalTag, jecTag, jecDB, muonDB = getEraConfiguration(era=options.era,isData=options.runOnData)

process = cms.Process("MiniAnalysis")

# Load the standard set of configuration modules
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')


# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag)


#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


# set input to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root'),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )
if options.runOnData:
    process.source.fileNames = cms.untracked.vstring('/store/data/Run2016B/MuonEG/MINIAOD/23Sep2016-v3/00000/024ADA16-1F98-E611-AD32-0242AC130005.root')

#apply lumi json, if passed in command line
if options.lumiJson:
    print 'Lumi sections will be selected with',options.lumiJson
    from FWCore.PythonUtilities.LumiList import LumiList
    myLumis = LumiList(filename = options.lumiJson).getCMSSWString().split(',')
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess.extend(myLumis)



#analysis
process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
process.analysis.muonCorrFile=cms.string(muonDB)
process.analysis.applyLeptonCorrections=cms.bool(options.applyLeptonCorrections)
process.analysis.saveTree=cms.bool(options.saveTree)
process.analysis.savePF=cms.bool(options.savePF)
print 'Ntuplizer configuration is as follows'
if process.analysis.applyLeptonCorrections:
    print '\t Muon corrections served from',process.analysis.muonCorrFile
else:
    print '\t Lepton corrections are disabled'
if not process.analysis.saveTree :
    print '\t Summary tree won\'t be saved'
if not process.analysis.savePF :
    print 'Summary PF info won\'t be saved'
    

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True)
)

#pseudo-top
if not options.runOnData:
    process.load('TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi')
    process.pseudoTop.leptonMinPt=cms.double(20)
    process.pseudoTop.leptonMaxEta=cms.double(2.5)
    process.pseudoTop.jetMaxEta=cms.double(5.0)

#EGM smearer https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer
process.load('Configuration.StandardSequences.Services_cff')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                                       engineName = cms.untracked.string('TRandom3'),
                                                                                       )
                                                   )
process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
process.calibratedPatElectrons.correctionFile = 'EgammaAnalysis/ElectronTools/data/ScalesSmearings/Winter_2016_reReco_v1_ele'
process.calibratedPatElectrons.isMC=False if options.runOnData else True
print 'Using calibrated electrons with corrections from',process.calibratedPatElectrons.correctionFile

# Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff'] 
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


#jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from TopLJets2015.TopAnalysis.customizeJetTools_cff import *
customizeJetTools(process=process,
                  jecTag=jecTag,
                  jecDB=jecDB,
                  baseJetCollection=options.baseJetCollection,
                  runOnData=options.runOnData)

#tfile service
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename)
                                   )

if options.runOnData:
    process.p = cms.Path(process.calibratedPatElectrons*process.egmGsfElectronIDSequence*process.analysis)
else:
    process.p = cms.Path(process.calibratedPatElectrons*process.egmGsfElectronIDSequence*process.pseudoTop*process.analysis)
