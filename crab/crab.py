########################
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8___coolv1'
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../runme_cfg.py'
config.Data.inputDBS = 'global'


config.Data.inputDataset =  '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 250000

config.Data.outLFNDirBase = '/store/user/theofil/test' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = ''
config.Site.storageSite = 'T2_CH_CERN'
