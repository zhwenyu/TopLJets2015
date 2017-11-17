from WMCore.Configuration import Configuration
import os
config = Configuration()

config.section_("General")
config.General.requestName = 'MC13TeV_TTJets'
config.General.workArea = 'my_grid'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/user/q/qhassan/workQhassan/TOPLJets2017/CMSSW_9_2_3/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py'
config.JobType.disableAutomaticOutputCollection = False
config.JobType.pyCfgParams = ['runOnData=False']
config.JobType.inputFiles = ['jec_MC.db']

config.section_("Data")
config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10_ext1-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
#config.Data.outLFNDirBase = '/store/user/TOP/' #% (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'November2017_MC_analysis'
config.Data.ignoreLocality = False

config.section_("Site")
config.Site.storageSite = 'T2_IT_Legnaro'

