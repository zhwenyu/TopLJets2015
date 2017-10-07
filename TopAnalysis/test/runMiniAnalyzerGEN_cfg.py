import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('outFilename', 'MiniEvents.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('inputFile', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input file to process"
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
options.register('nevents', 10, VarParsing.multiplicity.singleton, VarParsing.varType.int, "# events to generate")
options.register('scale', 1., VarParsing.multiplicity.singleton, VarParsing.varType.float, "factor for fsr ren scale")
options.register('asfsr', 0.1365, VarParsing.multiplicity.singleton, VarParsing.varType.float, "alpha_s for fsr")
options.register('asisr', 0.1180, VarParsing.multiplicity.singleton, VarParsing.varType.float, "alpha_s for isr")
options.register('me', 'on', VarParsing.multiplicity.singleton, VarParsing.varType.string, "ME corrections")
options.register('generator', 'pythia8', VarParsing.multiplicity.singleton, VarParsing.varType.string, "PS generator")
options.register('cr', 'default', VarParsing.multiplicity.singleton, VarParsing.varType.string, "color reconnection mode")
options.parseArguments()

process = cms.Process("MiniAnalysisGEN")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

# set input to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.nevents) )

process.options   = cms.untracked.PSet(
   # wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True)
)
# Input source
#process.source = cms.Source("PoolSource",
#    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
#    fileNames = cms.untracked.vstring( 
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E0E60B0A-EFD9-E411-944B-002590494C44.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E216B380-EFD9-E411-A286-003048FEC15C.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E41672D8-EFD9-E411-96A4-02163E00F2F9.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E41DFBBB-EFD9-E411-B7E3-02163E00F319.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E4F6FEEE-EED9-E411-ABA1-02163E00F4EF.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E612966B-EFD9-E411-A4DB-02163E00EAD0.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E85BB4DA-EED9-E411-8FE2-02163E00E9BC.root',
#    ),
#    inputCommands = cms.untracked.vstring('keep *', 
#        'drop LHEXMLStringProduct_*_*_*'),
#    secondaryFileNames = cms.untracked.vstring()
#)

process.source = cms.Source("EmptySource")
process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    nEvents = cms.untracked.uint32(options.nevents),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh'),
    numberOfParameters = cms.uint32(1),
    args = cms.vstring('/eos/user/m/mseidel/ReReco2016/b312177_merged/TT_hdamp_TuneT4_noweights_NNPDF30_13TeV_powheg_hvq.tgz')
)

if (options.generator == "pythia8"):
    from Configuration.Generator.Pythia8CommonSettings_cfi import *
    from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
    from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *
    process.generator = cms.EDFilter("Pythia8HadronizerFilter",
        maxEventsToPrint = cms.untracked.int32(0),
        pythiaPylistVerbosity = cms.untracked.int32(0),
        filterEfficiency = cms.untracked.double(1.0),
        pythiaHepMCVerbosity = cms.untracked.bool(False),
        comEnergy = cms.double(13000.),
        PythiaParameters = cms.PSet(
            pythia8CommonSettingsBlock,
            pythia8CUEP8M1SettingsBlock,
            pythia8PowhegEmissionVetoSettingsBlock,
            processParameters = cms.vstring(
                'POWHEG:nFinal = 2',
                'TimeShower:mMaxGamma = 1.0',
                'Tune:pp 14',
                'Tune:ee 7',
                'MultipartonInteractions:ecmPow=0.25208',
                'SpaceShower:alphaSvalue=0.1108',
                'PDF:pSet=LHAPDF6:NNPDF30_lo_as_0130',
                'MultipartonInteractions:pT0Ref=2.197139e+00',
                'MultipartonInteractions:expPow=1.608572e+00',
                'ColourReconnection:range=6.593269e+00',
            ),
            fsrParameters = cms.vstring(
                'TimeShower:renormMultFac   = ' + str(options.scale**2),
                'TimeShower:factorMultFac   = ' + str(options.scale**2),
                'TimeShower:alphaSvalue     = ' + str(options.asfsr),
                'TimeShower:MEcorrections   = ' + options.me,
            ),
            isrParameters = cms.vstring(
                'SpaceShower:alphaSvalue     = ' + str(options.asisr),
                ),
            parameterSets = cms.vstring('pythia8CommonSettings',
                #'pythia8CUEP8M1Settings',
                'pythia8PowhegEmissionVetoSettings',
                'processParameters',
                'fsrParameters'
            )
        )
    )

    if (options.cr == "off"):
        process.generator.PythiaParameters.processParameters.append('ColourReconnection:reconnect = off')

    # Tune from TOP-RunIISummer15wmLHEGS-00176
    if ('qcd' in options.cr):
        process.generator.PythiaParameters.processParameters = cms.vstring(
            'POWHEG:nFinal = 2', ## Number of final state particles
            ## (BEFORE THE DECAYS) in the LHE
            ## other than emitted extra parton
            'TimeShower:mMaxGamma = 1.0',#cutting off lepton-pair production
            ##in the electromagnetic shower
            ##to not overlap with ttZ/gamma* samples
            'Tune:pp 14',
            'Tune:ee 7',
            'MultipartonInteractions:ecmPow=0.25208',
            'SpaceShower:alphaSvalue=0.1108',
            'PDF:pSet=LHAPDF6:NNPDF30_lo_as_0130',
            #'MultipartonInteractions:pT0Ref=2.197139e+00',
            #'MultipartonInteractions:expPow=1.608572e+00',
            #'ColourReconnection:range=6.593269e+00',
            'ColourReconnection:mode = 1',
            'BeamRemnants:remnantMode=1',
            'StringZ:aLund = 0.38',
            'StringZ:bLund = 0.64',
            'StringFlav:probQQtoQ = 0.078',
            'StringFlav:probStoUD = 0.2',
            'StringFlav:probQQ1toQQ0join=0.0275,0.0275,0.0275',
            'ColourReconnection:allowDoubleJunRem = on', #TODO Not recommended?
            'ColourReconnection:timeDilationPar=15.86',
            'ColourReconnection:m0=1.204',
            'ColourReconnection:junctionCorrection=0.1222',
            'MultipartonInteractions:pT0Ref=2.174',
            'MultipartonInteractions:expPow=1.312',
        )
    
    # Tune from TOP-RunIISummer15wmLHEGS-00223
    if ('move' in options.cr):
        process.generator.PythiaParameters.processParameters = cms.vstring(
            'POWHEG:nFinal = 2', ## Number of final state particles
            ## (BEFORE THE DECAYS) in the LHE
            ## other than emitted extra parton
            'TimeShower:mMaxGamma = 1.0',#cutting off lepton-pair production
            ##in the electromagnetic shower
            ##to not overlap with ttZ/gamma* samples
            'Tune:pp 14',
            'Tune:ee 7',
            'MultipartonInteractions:ecmPow=0.25208',
            'SpaceShower:alphaSvalue=0.1108',
            'PDF:pSet=LHAPDF6:NNPDF30_lo_as_0130',
            'MultipartonInteractions:pT0Ref=2.30e+00',
            'MultipartonInteractions:expPow=1.346e+00',
            'ColourReconnection:mode=2',
            'ColourReconnection:m2Lambda=1.89',
          )
    
    if ('erd' in options.cr):
        process.generator.PythiaParameters.processParameters.append('PartonLevel:earlyResDec = on')


if (options.generator == "herwigpp"):
    from Configuration.Generator.HerwigppDefaults_cfi import *
    from Configuration.Generator.HerwigppUE_EE_5C_cfi import *
    from Configuration.Generator.HerwigppPDF_CTEQ6_LO_cfi import * # Import CTEQ6L PDF as shower pdf
    from Configuration.Generator.HerwigppPDF_NNPDF30_NLO_cfi import herwigppPDFSettingsBlock as herwigppHardPDFSettingsBlock 	# Import NNPDF30 NLO as PDF of the hard subprocess
    from Configuration.Generator.HerwigppEnergy_13TeV_cfi import *
    from Configuration.Generator.HerwigppLHEFile_cfi import *
    from Configuration.Generator.HerwigppMECorrections_cfi import *

    process.generator = cms.EDFilter("ThePEGHadronizerFilter",
	      herwigDefaultsBlock,
	      herwigppUESettingsBlock,
	      herwigppPDFSettingsBlock,
	      herwigppHardPDFSettingsBlock,			# Implementing renamed NNPDF30 config block
	      herwigppEnergySettingsBlock,
	      herwigppLHEFileSettingsBlock,
	      herwigppMECorrectionsSettingsBlock,

	      configFiles = cms.vstring(),
	      parameterSets = cms.vstring(
		        'hwpp_cmsDefaults',
		        'hwpp_ue_EE5C',
		        'hwpp_cm_13TeV',
		        'hwpp_pdf_CTEQ6L1',			# Shower PDF matching with the tune
		        'hwpp_pdf_NNPDF30NLO_Hard',		# PDF of hard subprocess
		        'hwpp_LHE_Powheg_DifferentPDFs',
		        'hwpp_MECorr_' + options.me.capitalize(),
		        'productionParameters',
	      ),
	      productionParameters = cms.vstring(
		        'set /Herwig/Shower/AlphaQCD:RenormalizationScaleFactor ' + str(options.scale),
		        'set /Herwig/Shower/AlphaQCD:AlphaMZ ' + str(options.asfsr)
        ),
        crossSection = cms.untracked.double(-1),
        filterEfficiency = cms.untracked.double(1.0),
        #dumpConfig  = cms.untracked.string("dump_me_" + options.me + ".config"),
    )

#pseudo-top
process.load("TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi")
process.pseudoTop.src = 'generator:unsmeared'
process.pseudoTop.runTopReconstruction = False

process.load('TopQuarkAnalysis.BFragmentationAnalyzer.bfragWgtProducer_cfi')

#tfile service
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename)
                                   )


#analysis
#process.analysis = cms.EDAnalyzer("MiniAnalyzerGEN",
#                          saveTree               = cms.bool(True),
#                          savePF                 = cms.bool(True),
#                          )

process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
process.analysis.saveTree = options.saveTree
process.analysis.savePF   = options.savePF
process.analysis.runOnGEN = True
process.analysis.prunedGenParticles = 'genParticles'
if not process.analysis.saveTree :
    print '\t Summary tree won\'t be saved'
if not process.analysis.savePF :
    print 'Summary PF info won\'t be saved'

# Path and EndPath definitions
#process.lhe_step = cms.Path(process.externalLHEProducer)
process.generation_step = cms.Path(process.generator*process.pseudoTop*process.bfragWgtProducer*process.analysis)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# FastGenParticleCandidateProducer broken, replace by GenParticleProducer
# (actually no idea where this is from)
process.genParticleCandidates = cms.EDProducer("GenParticleProducer",
    abortOnUnknownPDGCode = cms.untracked.bool(False),
    saveBarCodes = cms.untracked.bool(True),
    src = cms.InputTag("generatorSmeared")
)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step)
# filter all path with the production filter sequence
for path in process.paths:
	  if path in ['lhe_step']: continue
	  getattr(process,path)._seq = process.generator * getattr(process,path)._seq 
