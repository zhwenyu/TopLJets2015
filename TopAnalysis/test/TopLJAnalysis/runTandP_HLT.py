import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('func',
		 'BWResCBExp',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "function to fit"
                 )
options.register('tree',
		 'pPb_8.16TeV',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "tree to use"
                 )
options.register('data',
		 True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "is data"
                 )
options.parseArguments()

FITFUNC = options.func
OUTFILE='hltfits_%s_%d_%s.root'%(options.tree,options.data,options.func)
from tandpTrees import FILES
INFILE=FILES[options.tree][0] if options.data else FILES[options.tree][1]

print FITFUNC,INFILE,OUTFILE

process.TnP_Muon_ID = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    ## Input, output 
    InputFileNames = cms.vstring(INFILE),
    OutputFileName = cms.string(OUTFILE),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    #WeightVariable = cms.string("weight"),
    ## Variables for binning
    Variables = cms.PSet(
        mass   = cms.vstring("Tag-Probe Mass", "60", "120", "GeV"),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV"),
        eta    = cms.vstring("Probe #eta", "-2.1", "2.1", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        combRelIsoPF04dBeta = cms.vstring("PF Combined Relative Iso", "-100", "99999", "")
    ),
                                     
    ## Flags you want to use to define numerator and possibly denominator
   Categories = cms.PSet(
	tag_hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12 = cms.vstring("tag_hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12", "dummy[pass=1,fail=0]"),
        hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12     = cms.vstring("hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12", "dummy[pass=1,fail=0]"),
	Tight2012                                         = cms.vstring("Tight2012", "dummy[pass=1,fail=0]"),
    ),

                                     
    ## What to fit
    Efficiencies = cms.PSet(
        MuID_pt = cms.PSet(
            UnbinnedVariables = cms.vstring("mass"),
            EfficiencyCategoryAndState = cms.vstring("hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12", "pass"),
            BinnedVariables = cms.PSet(
                ## Binning in continuous variables
                eta = cms.vdouble(-2.1, 2.1),
                pt = cms.vdouble( 20, 30, 40, 50, 100 ),
                ## flags and conditions required at the denominator, 
                tag_hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12 = cms.vstring("pass"),
                Tight2012 = cms.vstring('pass'),
                combRelIsoPF04dBeta = cms.vdouble(0,0.15)
            ),
            BinToPDFmap = cms.vstring(FITFUNC), ## PDF to use, as defined below
        ),
        MuID_eta = cms.PSet(
            UnbinnedVariables = cms.vstring("mass"),
            EfficiencyCategoryAndState = cms.vstring("hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12", "pass"),
            BinnedVariables = cms.PSet(
                abseta = cms.vdouble(0, 0.6, 1.1, 1.6, 2.1),
                pt = cms.vdouble(30, 100),
                tag_hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12 = cms.vstring("pass"),
                Tight2012 = cms.vstring('pass'),
                combRelIsoPF04dBeta = cms.vdouble(0,0.15)
            ),
            BinToPDFmap = cms.vstring(FITFUNC), ## PDF to use, as defined below
        ),
    ),

    ## PDF for signal and background (double voigtian + exponential background)
    PDFs = cms.PSet(
        CB2Exp = cms.vstring(
            "CBShape::signalPass(mass, peakPass[90,80,100], sigmaPass[2.5,1,20], alphaPass[0.5,0,3], nPass[5,1,8])",
            "CBShape::signalFail(mass, peakFail[90,80,100], sigmaFail[2.5,1,20], alphaFail[0.5,0,3], nFail[5,1,8])",
            "RooExponential::backgroundPass(mass, lp[0,-5,5])",
            "RooExponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9,0.5,1.0]"
            ),
	VoigtExp = cms.vstring(
		"Voigtian::signal(mass, mean[91,85,95], width[3,1,10], sigma[3,1,10])",
		"Exponential::backgroundPass(mass, lp[0,-5,5])",
		"Exponential::backgroundFail(mass, lf[0,-5,5])",
		"efficiency[0.9,0,1]",
		"signalFractionInPassing[0.9]"
	),
	BWResCBExp = cms.vstring(
            "BreitWigner::bw(mass, m0[91.2,81.2,101.2], width[2.495,1,10])",	
            "RooCBShape::res(mass, peak[0], sigma[1.7,1,10], alpha[1.8,0,3], n[0.8,0,10])",
            "FCONV::signal(mass, bw, res)",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0.5,1]",
            "signalFractionInPassing[0.9]",
            ),
	BWResCBCheb = cms.vstring(
            "BreitWigner::bw(mass, m0[91.2,81.2,101.2], width[2.495,1,10])",            
            "RooCBShape::res(mass, peak[0], sigma[1.7,1,10], alpha[1.8,0,3], n[0.8,0,10])",
            "FCONV::signal(mass, bw, res)",
            "Chebychev::backgroundPass(mass, {c1p[0,-10,10], c2p[0,-10,10], c3p[0,-10,10]})",
            "Chebychev::backgroundFail(mass, {c1f[0,-10,10], c2f[0,-10,10], c3f[0,-10,10]})",
            "efficiency[0.9,0.5,1]",
            "signalFractionInPassing[0.9,0.7,1.0]",
            ),
        ),

    ## How to do the fit
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),
    saveDistributionsPlot = cms.bool(True),
    NumCPU = cms.uint32(1), ## leave to 1 for now, RooFit gives funny results otherwise
    SaveWorkspace = cms.bool(True),
)

process.p = cms.Path(process.TnP_Muon_ID)
