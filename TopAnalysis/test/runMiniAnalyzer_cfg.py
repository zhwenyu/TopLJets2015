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
options.register('era', 'era2018',
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

process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi") 


#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#EGM
from TopLJets2015.TopAnalysis.customizeEGM_cff import customizeEGM
try:
      customizeEGM(process=process,runOnData=options.runOnData)
except:
      print 'Failed to customize EGM sequence, will use default in MINIAOD'
      options.useRawLeptons=True
      
# set input to process
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIISummer17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v3/10000/00F1AD27-CF99-E711-A2F0-0CC47AC087AE.root'),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )

if options.runOnData:
      process.source.fileNames = cms.untracked.vstring('/store/data/Run2018A/SingleMuon/MINIAOD/PromptReco-v1/000/315/252/00000/A4345B88-404B-E811-894C-FA163EECC2C3.root')

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

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                    inputPruned = cms.InputTag("prunedGenParticles"),
                        inputPacked = cms.InputTag("packedGenParticles"),
)
from GeneratorInterface.RivetInterface.genParticles2HepMC_cfi import genParticles2HepMC
process.genParticles2HepMC = genParticles2HepMC.clone( genParticles = cms.InputTag("mergedGenParticles") )
#process.load('TopQuarkAnalysis.BFragmentationAnalyzer.bfragWgtProducer_cfi')

#apply lumi json, if passed in command line
if options.lumiJson:
    print 'Lumi sections will be selected with',options.lumiJson
    from FWCore.PythonUtilities.LumiList import LumiList
    myLumis = LumiList(filename = options.lumiJson).getCMSSWString().split(',')
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess.extend(myLumis)


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

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename)
                                  )

process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
from TopLJets2015.TopAnalysis.miniAnalyzer_cfi import triggerLists
print 'Ntuplizer configuration is as follows'
process.analysis.saveTree=cms.bool(options.saveTree)
process.analysis.savePF=cms.bool(options.savePF)
process.analysis.useRawLeptons=cms.bool(options.useRawLeptons)
if options.useRawLeptons:
      print 'Switched off corrections for leptons'
      try:
            process.selectedElectrons.src=cms.InputTag('slimmedElectrons')
            process.analysis.electrons=cms.InputTag('selectedElectrons')
      except:
            process.analysis.electrons=cms.InputTag('slimmedElectrons')
            
if not process.analysis.saveTree :
    print '\t Summary tree won\'t be saved'
if not process.analysis.savePF :
    print 'Summary PF info won\'t be saved'
eraNumber= int(filter(str.isdigit, options.era))
process.analysis.triggersToUse=triggerLists[eraNumber]
print 'Applying %d trigger list'%eraNumber

if options.runOnData:
      process.analysis.metFilterBits = cms.InputTag("TriggerResults","","RECO")

process.ana_step=cms.Path(process.analysis)

#schedule
try:
      process.custom_step=cms.Path(process.egammaScaleSmearAndVIDSeq)
      if options.runOnData:
            process.schedule=cms.Schedule(process.custom_step,process.ana_step)
      else:
            process.gen_step=cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel)
            process.schedule=cms.Schedule(process.custom_step,process.gen_step,process.ana_step)
except:
      print 'Failed to setup standard scheduler, will just read what is available in MINIAOD'
      if options.runOnData:
            process.schedule=cms.Schedule(process.ana_step)
      else:
            process.gen_step=cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel)
            process.schedule=cms.Schedule(process.gen_step,process.ana_step)
