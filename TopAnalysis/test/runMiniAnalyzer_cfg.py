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

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
#                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff']
                 

#add them to the VID producer
for idmod in my_id_modules:
      setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# set input to process
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIISummer17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v3/10000/00F1AD27-CF99-E711-A2F0-0CC47AC087AE.root',
                            '/store/mc/RunIISummer17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v3/10000/04582658-7A99-E711-8C87-90B11C2AA430.root',
                            '/store/mc/RunIISummer17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v3/10000/06320BFA-D099-E711-9E2B-001EC9AE27F4.root',
                            '/store/mc/RunIISummer17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v3/10000/105E5310-9799-E711-8BED-0025905A60D2.root',
                            '/store/mc/RunIISummer17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v3/10000/18ABEEC9-CE99-E711-97F0-0CC47A6C1866.root',
                            '/store/mc/RunIISummer17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v3/10000/1AE49EAA-D099-E711-9E94-48D539D33363.root',
                            '/store/mc/RunIISummer17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v3/10000/1C2D020E-9799-E711-97CE-0025905C53B2.root',
                            '/store/mc/RunIISummer17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v3/10000/2A0CA3BC-8299-E711-86EF-0242AC11000E.root',
                            '/store/mc/RunIISummer17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v3/10000/2C373EF3-CE99-E711-B09F-90E2BAC9B7A8.root',
                            '/store/mc/RunIISummer17MiniAOD/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/92X_upgrade2017_realistic_v10-v3/10000/36892FF1-CE99-E711-8539-7CD30ACE1732.root'),
                            
                            #fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root'
#                                                              ),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )

if options.runOnData:
      process.source.fileNames = cms.untracked.vstring('/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/678/00000/7210909D-AA66-E711-982B-02163E0146EB.root',
                                                       '/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/679/00000/985E2CD1-A766-E711-BA46-02163E01421A.root',
                                                       '/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/681/00000/84C1155B-AE66-E711-9C1C-02163E0134E0.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/853/00000/245298CA-7F68-E711-BCB5-02163E01A505.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/853/00000/529C82FB-8768-E711-B2EB-02163E0142F7.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/853/00000/60253CBC-A368-E711-A664-02163E0144F7.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/853/00000/7E4FC127-8468-E711-AFE9-02163E0144D8.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/853/00000/F0A2A315-8F68-E711-859E-02163E019E36.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/854/00000/B6958769-6768-E711-A751-02163E011A09.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/855/00000/98FD21AF-6F68-E711-BB0C-02163E013390.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/996/00000/12B43734-166A-E711-81E5-02163E011C1F.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/996/00000/3A9296E0-186A-E711-BC21-02163E01A3A5.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/996/00000/6C27785A-146A-E711-9C7C-02163E01A638.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/996/00000/A2C3BC09-356A-E711-A65C-02163E0135EF.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/996/00000/D0F7F995-226A-E711-BE9B-02163E01A382.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/996/00000/EE20CC73-1C6A-E711-B058-02163E019E36.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/997/00000/2650595E-2E6A-E711-8A6C-02163E011D43.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/298/998/00000/469068AD-156A-E711-B6EE-02163E01A505.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/000/00000/AC28272D-1F6A-E711-83D6-02163E013490.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/000/00000/C66DA99F-256A-E711-90C2-02163E0137FC.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/042/00000/3C0EBDA6-8A6A-E711-BD36-02163E01A5D4.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/0C47A1A3-B96A-E711-B9C0-02163E019C3E.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/1E53C8AC-AC6A-E711-96C3-02163E014235.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/203FA7D4-C16A-E711-B7EC-02163E01479A.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/5CBB654C-AB6A-E711-9A73-02163E0135F2.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/5E3E5850-AD6A-E711-8400-02163E011E56.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/6C89BB14-AF6A-E711-959C-02163E011926.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/72B4F23E-BD6A-E711-A7CA-02163E011EE8.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/7A646BBE-B26A-E711-9606-02163E011C1F.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/88089B90-CC6A-E711-9D68-02163E0137FC.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/983B0115-B86A-E711-8BB2-02163E0133AF.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/9ACEF047-B06A-E711-BE99-02163E0141FB.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/A0C79BF0-AD6A-E711-A93A-02163E01A619.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/CAE0F63C-B46A-E711-A180-02163E01A208.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/D2D8092D-B66A-E711-8E17-02163E01A6A3.root',
'/store/data/Run2017B/SingleElectron/MINIAOD/PromptReco-v2/000/299/061/00000/FC00A170-B16A-E711-B9D5-02163E0129A0.root')                                  
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                    inputPruned = cms.InputTag("prunedGenParticles"),
                        inputPacked = cms.InputTag("packedGenParticles"),
)
from GeneratorInterface.RivetInterface.genParticles2HepMC_cfi import genParticles2HepMC
process.genParticles2HepMC = genParticles2HepMC.clone( genParticles = cms.InputTag("mergedGenParticles") )
#process.load("TopQuarkAnalysis.TopEventProducers.producers.particleLevel_cfi")
#process.load('TopQuarkAnalysis.BFragmentationAnalyzer.bfragWgtProducer_cfi')

#apply lumi json, if passed in command line
if options.lumiJson:
    print 'Lumi sections will be selected with',options.lumiJson
    from FWCore.PythonUtilities.LumiList import LumiList
    myLumis = LumiList(filename = options.lumiJson).getCMSSWString().split(',')
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess.extend(myLumis)

#EGM
#from TopLJets2015.TopAnalysis.customizeEGM_cff import *
#customizeEGM(process=process,runOnData=options.runOnData)
#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#                                                   calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
#                                                                                      engineName = cms.untracked.string('TRandom3'),
#                                                                                     )
#                                                   )
'''
process.source = cms.Source("EmptySource")
process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    nEvents = cms.untracked.uint32(options.nevents),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh'),
    numberOfParameters = cms.uint32(1),
    args = cms.vstring('/eos/user/m/mseidel/ReReco2016/b312177_merged/TT_hdamp_TuneT4_noweights_NNPDF30_13TeV_powheg_hvq.tgz')
)
'''

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
                                   fileName = cms.string('histo.root')
                                  )

process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
print 'Ntuplizer configuration is as follows'
process.analysis.saveTree=cms.bool(options.saveTree)
process.analysis.savePF=cms.bool(options.savePF)
if options.useRawLeptons:
    process.selectedElectrons.src=cms.InputTag('slimmedElectrons')
    process.analysis.electrons=cms.InputTag('selectedElectrons')
    process.analysis.useRawLeptons=cms.bool(True)
    print 'Switched off corrections for leptons'
if not process.analysis.saveTree :
    print '\t Summary tree won\'t be saved'
if not process.analysis.savePF :
    print 'Summary PF info won\'t be saved'

if options.runOnData:
    process.analysis.metFilterBits = cms.InputTag("TriggerResults","","RECO")

if options.runOnData:
    process.p = cms.Path(process.analysis*process.egmGsfElectronIDSequence)
else:
    process.p = cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel*process.egmGsfElectronIDSequence*process.analysis)

