import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('standard')
options.register('yodafile',  'test.yoda', VarParsing.multiplicity.singleton, VarParsing.varType.string, "Name of yoda output file")
options.register('scale',     1.,          VarParsing.multiplicity.singleton, VarParsing.varType.float,  "factor for fsr ren scale")
options.register('asfsr',     0.1365,      VarParsing.multiplicity.singleton, VarParsing.varType.float,  "alpha_s for fsr")
options.register('asisr',     0.1180,      VarParsing.multiplicity.singleton, VarParsing.varType.float,  "alpha_s for isr")
options.register('me',        'on',        VarParsing.multiplicity.singleton, VarParsing.varType.string, "ME corrections")
options.register('generator', 'pythia8',   VarParsing.multiplicity.singleton, VarParsing.varType.string, "PS generator")
options.register('cr',        'default',   VarParsing.multiplicity.singleton, VarParsing.varType.string, "color reconnection mode")
if(hasattr(sys, "argv")):
    options.parseArguments()
print options

process = cms.Process("runRivetAnalysis")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

process.source = cms.Source("EmptySource")

process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
                                             args = cms.vstring('/eos/user/m/mseidel/ReReco2016/b312177_merged/TT_hdamp_TuneT4_noweights_NNPDF30_13TeV_powheg_hvq.tgz'),
                                             #args = cms.vstring('/cvmfs/cms.cern.ch/phys_generator/gridpacks/slc6_amd64_gcc481/13TeV/powheg/V2/TT_hdamp_Tune4_NNPDF30/TT_hdamp_Tune4_NNPDF30_13TeV_powheg_hvq.tgz'),
                                             nEvents = cms.untracked.uint32(options.maxEvents),
                                             numberOfParameters = cms.uint32(1),
                                             outputFile = cms.string('cmsgrid_final.lhe'),
                                             scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
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

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2017_TOP_17_015')
process.rivetAnalyzer.OutputFile      = options.yodafile
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection    = 831.76 # NNLO (arXiv:1303.6254)

#pseudo-top
process.load("TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi")
process.pseudoTop.src = 'generator:unsmeared'
process.pseudoTop.runTopReconstruction = False

process.load('TopQuarkAnalysis.BFragmentationAnalyzer.bfragWgtProducer_cfi')

#tfile service
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.yodafile.replace('.yoda','.root'))
                                   )

process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
process.analysis.skim     = False
process.analysis.saveTree = True
process.analysis.savePF   = True
process.analysis.runOnGEN = True
process.analysis.prunedGenParticles = 'genParticles'

process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.genParticles.src='generator:unsmeared'

process.p = cms.Path(process.externalLHEProducer*process.generator*process.genParticles*
                     process.rivetAnalyzer*
                     process.pseudoTop*process.bfragWgtProducer*process.analysis)

