import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnData', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run this on real data"
                 )
options.register('noParticleLevel', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Do not run the particleLevel sequence"
                 )
options.register('era', 'era2017',
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
options.register('savePF', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'save PF candidates'
                 )
options.parseArguments()

#start process
process = cms.Process("MiniAnalysis")

#get the configuration to apply
from TopLJets2015.TopAnalysis.EraConfig import getEraConfiguration
globalTag, jecTag, jecDB, jerTag, jerDB = getEraConfiguration(era=options.era,isData=options.runOnData)

# Load the standard set of configuration modules
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#EGM customization
from TopLJets2015.TopAnalysis.customizeEGM_cff import customizeEGM
process.egammaPostRecoSeq=customizeEGM(process=process,era=options.era)

# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag)

process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi") 

# particle level definitions
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                                            inputPruned = cms.InputTag("prunedGenParticles"),
                                            inputPacked = cms.InputTag("packedGenParticles"),
                                            )
process.load('GeneratorInterface.RivetInterface.genParticles2HepMC_cfi')
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")

#jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from TopLJets2015.TopAnalysis.customizeJetTools_cff import *
customizeJetTools(process=process,
                  jecDB=jecDB,
                  jecTag=jecTag,
                  jerDB=jerDB,
                  jerTag=jerTag,
                  baseJetCollection=options.baseJetCollection,
                  runOnData=options.runOnData)

#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# set input to process
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FE7ABAEB-4A42-E811-87A3-0CC47AD98D26.root'),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )
if options.runOnData:
      process.source.fileNames = cms.untracked.vstring('')

if options.inputFile:
      fileList=[]
      if '.root' in options.inputFile :
            fileList=[options.inputFile]
      else:
            import os
            inDir=options.inputFile
            fileList = ['file:'+os.path.join(inDir,f) for f in os.listdir(inDir)]
      process.source.fileNames = cms.untracked.vstring(fileList)
print  "Processing",process.source.fileNames

#apply lumi json, if passed in command line
if options.lumiJson:
    print 'Lumi sections will be selected with',options.lumiJson
    from FWCore.PythonUtilities.LumiList import LumiList
    myLumis = LumiList(filename = options.lumiJson).getCMSSWString().split(',')
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess.extend(myLumis)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename)
                                  )

#analysis
process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
process.analysis.saveTree=cms.bool(options.saveTree)
process.analysis.savePF=cms.bool(options.savePF)
process.analysis.jetIdToUse=cms.string('tightLepVeto' if 'era2017' in options.era else 'looseID')
if not process.analysis.saveTree :
    print '\t Summary tree won\'t be saved'
if not process.analysis.savePF :
    print 'Summary PF info won\'t be saved'

if options.runOnData:
      process.analysis.metFilterBits = cms.InputTag("TriggerResults","","RECO")
#process.ana_step=cms.Path(process.analysis)

#schedule execution
toSchedule=[]
if process.egammaPostRecoSeq:
      process.egm=cms.Path(process.egammaPostRecoSeq) 
      toSchedule.append( process.egm )
if not (options.runOnData or options.noParticleLevel):
      process.mctruth=cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel)
      toSchedule.append( process.mctruth )
process.ana=cms.Path(process.analysis)
toSchedule.append( process.ana )

process.schedule=cms.Schedule( (p for p in toSchedule) )


#if options.runOnData or options.noParticleLevel:
#      if process.egammaPostRecoSeq:
#            process.schedule=cms.Schedule(cms.Path(process.egammaPostRecoSeq),process.ana_step)
#      else:
#            process.schedule=cms.Schedule(process.ana_step)
#else:
#      process.gen_step=cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel)
#      if process.egammaPostRecoSeq:
#            process.schedule=cms.Schedule(cms.Path(process.egammaPostRecoSeq),process.gen_step,process.ana_step)
#      else:
#            process.schedule=cms.Schedule(process.gen_step,process.ana_step)
