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
options.register('baseJetCollection','slimmedJets',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Base jet collection"
                 )
options.register('inputFile', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input file to process"
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
options.parseArguments()

#get the configuration to apply
from TopLJets2015.TopAnalysis.EraConfig import getEraConfiguration
globalTag, jecTag, jecDB = getEraConfiguration(era=options.era,isData=options.runOnData)

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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root'),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )
if options.runOnData:
    process.source.fileNames = cms.untracked.vstring('/store/data/Run2016B/MuonEG/MINIAOD/23Sep2016-v3/120000/E8C2C234-BF97-E611-82BB-0025907B4EE4.root')

if options.inputFile and len(options.inputFile)!=0:
    process.source.fileNames = cms.untracked.vstring(options.inputFile)

#apply lumi json, if passed in command line
if options.lumiJson:
    print 'Lumi sections will be selected with',options.lumiJson
    from FWCore.PythonUtilities.LumiList import LumiList
    myLumis = LumiList(filename = options.lumiJson).getCMSSWString().split(',')
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess.extend(myLumis)

process.options   = cms.untracked.PSet(
   # wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True)
)

#pseudo-top
if not options.runOnData:
    process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
    process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
        inputPruned = cms.InputTag("prunedGenParticles"),
        inputPacked = cms.InputTag("packedGenParticles"),
    )

    from GeneratorInterface.RivetInterface.genParticles2HepMC_cfi import genParticles2HepMC
    process.genParticles2HepMC = genParticles2HepMC.clone( genParticles = cms.InputTag("mergedGenParticles") )

    process.load("TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi")
    process.pseudoTop.leptonMinPt=cms.double(20)
    process.pseudoTop.leptonMaxEta=cms.double(2.5)
    process.pseudoTop.jetMaxEta=cms.double(5.0)

    process.load('TopQuarkAnalysis.BFragmentationAnalyzer.bfragWgtProducer_cfi')

#EGM
from TopLJets2015.TopAnalysis.customizeEGM_cff import *
customizeEGM(process=process,runOnData=options.runOnData)
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                                       engineName = cms.untracked.string('TRandom3'),
                                                                                       )
                                                   )

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


#analysis
process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
process.analysis.saveTree=cms.bool(options.saveTree)
process.analysis.savePF=cms.bool(options.savePF)
process.analysis.egmscalecorr = process.calibratedPatElectrons.correctionFile
print 'Ntuplizer configuration is as follows'
if not process.analysis.saveTree :
    print '\t Summary tree won\'t be saved'
if not process.analysis.savePF :
    print 'Summary PF info won\'t be saved'
if options.runOnData:
    process.analysis.metFilterBits = cms.InputTag("TriggerResults","","RECO")


if options.runOnData:
    process.p = cms.Path(process.analysis)
    #process.p = cms.Path(process.egmSeq*process.jetmetSeq*process.analysis)
else:
    process.p = cms.Path(process.pseudoTop*process.bfragWgtProducer*process.analysis)
    #process.p = cms.Path(process.egmSeq*process.jetmetSeq*process.pseudoTop*process.analysis)
