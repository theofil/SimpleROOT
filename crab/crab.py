########################
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ZZTo4L_13TeV_powheg_pythia8___RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1___v1d'
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../runme_cfg.py'
config.Data.inputDBS = 'global'


config.Data.inputDataset =  '/ZZTo4L_13TeV_powheg_pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 250000

config.Data.outLFNDirBase = '/store/user/theofil/test' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = ''
config.Site.storageSite = 'T2_CH_CERN'
