########################
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'MuonEG_v1c'
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../runme_cfg.py'
config.Data.inputDBS = 'global'

config.Data.inputDataset = '/MuonEG/Run2015B-PromptReco-v1/MINIAOD'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 250000
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251642_13TeV_PromptReco_Collisions15_JSON.txt'

config.Data.outLFNDirBase = '/store/user/theofil/test' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = ''
config.Site.storageSite = 'T2_CH_CERN'
