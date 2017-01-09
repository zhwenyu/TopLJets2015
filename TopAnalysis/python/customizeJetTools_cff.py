import FWCore.ParameterSet.Config as cms

def customizeJetTools(process,jecTag,baseJetCollection,runOnData):

	#general configurations
	jecTag += '_DATA' if runOnData else '_MC'
	payload='AK4PFchs'
	jecLevels=['L1FastJet','L2Relative','L3Absolute']
	if 'Puppi' in baseJetCollection: 
		payload='AK4PFPuppi'
		jecLevels = ['L2Relative', 'L3Absolute']		
	if runOnData : jecLevels.append( 'L2L3Residual' )
	print '[customizeJetTools]',jecTag,payload,jecLevels

	#setup the source for JEC 
	process.load('CondCore.CondDB.CondDB_cfi')
	from CondCore.CondDB.CondDB_cfi import CondDB
	process.jec = cms.ESSource("PoolDBESSource",
				   DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
				   timetype = cms.string('runnumber'),
				   toGet = cms.VPSet( cms.PSet(record = cms.string('JetCorrectionsRecord'),
							       tag    = cms.string('JetCorrectorParametersCollection_%s_AK4PFchs'%jecTag),
							       label  = cms.untracked.string('AK4PFchs')
							       ),
						      cms.PSet(record = cms.string('JetCorrectionsRecord'),
							       tag    = cms.string('JetCorrectorParametersCollection_%s_AK4PFPuppi'%jecTag),
							       label  = cms.untracked.string('AK4PFPuppi')
							       ),
						      ## here you add as many jet types as you need
						      ## note that the tag name is specific for the particular sqlite file 
						      ), 
				   connect = cms.string('sqlite_file:data/era2016/%s.db'%jecTag)
				   )

	## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
	process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')


	
	from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
	updateJetCollection(
		process,
		labelName='UpdatedJECBTag',
		jetSource = cms.InputTag(baseJetCollection),
		jetCorrections = (payload, cms.vstring(jecLevels), 'None'),
		btagDiscriminators = ['deepFlavourJetTags:probudsg', 
				      'deepFlavourJetTags:probb', 
				      'deepFlavourJetTags:probc', 
				      'deepFlavourJetTags:probbb',
				      'deepFlavourJetTags:probcc',
				      'pfSimpleSecondaryVertexHighEffBJetTags',
				      'pfCombinedInclusiveSecondaryVertexV2BJetTags'
				      ],
		btagInfos = ['pfInclusiveSecondaryVertexFinderTagInfos']
		)

	process.updatedPatJetsSeq = cms.Sequence(process.patJetCorrFactorsUpdatedJECBTag*process.updatedPatJetsUpdatedJECBTag)
