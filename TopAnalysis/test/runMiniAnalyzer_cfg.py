import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnData', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run this on real data"
                 )
options.register('useRawLeptons', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Do not correct electrons/muons with smearing/energy scales"
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

process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi") 


#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000


# set input to process
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/PhaseISpring17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/FlatPU28to62_90X_upgrade2017_realistic_v20-v1/00000/0260246D-3129-E711-8B3A-001E67DDD0AA.root',
                            '/store/mc/PhaseISpring17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/FlatPU28to62_90X_upgrade2017_realistic_v20-v1/00000/02781287-E22A-E711-8EF8-A0000420FE80.root',
                            '/store/mc/PhaseISpring17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/FlatPU28to62_90X_upgrade2017_realistic_v20-v1/00000/0405C65F-3229-E711-9FEE-70106F4A469C.root'),
                            
                            #fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root'),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )

if options.runOnData:
      process.source.fileNames = cms.untracked.vstring('/store/data/Run2017A/SingleMuon/MINIAOD/PromptReco-v2/000/296/168/00000/084C505D-784C-E711-8140-02163E019DA4.root',
      '/store/data/Run2017A/SingleMuon/MINIAOD/PromptReco-v2/000/296/172/00000/66D29210-674C-E711-AFD2-02163E01A270.root',
      '/store/data/Run2017A/SingleMuon/MINIAOD/PromptReco-v2/000/296/173/00000/4C85F093-654C-E711-99AA-02163E01A203.root')
                                  
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                    inputPruned = cms.InputTag("prunedGenParticles"),
                        inputPacked = cms.InputTag("packedGenParticles"),
)
from GeneratorInterface.RivetInterface.genParticles2HepMC_cfi import genParticles2HepMC
process.genParticles2HepMC = genParticles2HepMC.clone( genParticles = cms.InputTag("mergedGenParticles") )
#process.load("TopQuarkAnalysis.TopEventProducers.producers.particleLevel_cfi")
#process.load('TopQuarkAnalysis.BFragmentationAnalyzer.bfragWgtProducer_cfi')

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histo.root')
                                  )

process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')

if options.runOnData:
    process.analysis.metFilterBits = cms.InputTag("TriggerResults","","RECO")

if options.runOnData:
    process.p = cms.Path(process.analysis)
else:
    process.p = cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel*process.analysis)

