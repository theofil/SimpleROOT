########################
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'SingleMuon___Run2015D-PromptReco-v4___coolv1'
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../runme_cfg.py'
config.Data.inputDBS = 'global'


config.Data.inputDataset =  '/SingleMuon/Run2015D-PromptReco-v4/MINIAOD'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 250000

config.Data.outLFNDirBase = '/store/user/theofil/test' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = ''
config.Site.storageSite = 'T2_CH_CERN'
