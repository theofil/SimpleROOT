########################
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ttbar_synch_v4_withISO'
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../runme_cfg.py'

config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset  = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.totalUnits = 25
config.Data.outLFN = '/store/user/theofil/test' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = ''

config.Site.storageSite = 'T2_CH_CERN'

config.JobType.allowNonProductionCMSSW = True ##  CMSSW_7_3_2 on slc6_amd64_gcc481 is not among supported releases 
