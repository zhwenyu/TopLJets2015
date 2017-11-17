from WMCore.Configuration import Configuration
import os
config = Configuration()

config.section_("General")
config.General.requestName = "Data13TeV_DoubleMuon_2017BPromptv2"
config.General.workArea = "my_grid"
config.General.transferOutputs=True

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "/afs/cern.ch/user/q/qhassan/workQhassan/TOPLJets2017/CMSSW_9_2_3/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py"
config.JobType.disableAutomaticOutputCollection = False
config.JobType.pyCfgParams = ['runOnData=True']
config.JobType.inputFiles = ['jec_DATA.db']

config.section_("Data")
config.Data.inputDataset = "/DoubleMuon/Run2017B-PromptReco-v2/MINIAOD"
config.Data.inputDBS = "global"
config.Data.splitting = "LumiBased"
config.Data.unitsPerJob = 5
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-302343_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.publication = True
config.Data.ignoreLocality = False
#config.Data.outLFNDirBase = 'None/6c5beb6/'

config.section_("Site")
config.Site.storageSite = "T2_IT_Legnaro"
