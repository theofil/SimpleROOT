########################
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'DoubleMuon_v1'
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../runme_cfg.py'

#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'
#config.Data.inputDataset = '/DoubleEG/Run2015B-PromptReco-v1/MINIAOD'
config.Data.inputDataset = '/DoubleMuon/Run2015B-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/MuonEG/Run2015B-PromptReco-v1/MINIAOD'

#
#   /EGamma/Run2015B-PromptReco-v1/MINIAOD
#   /MET/Run2015B-PromptReco-v1/MINIAOD
#   /SingleElectron/Run2015B-PromptReco-v1/MINIAOD
#   /SingleMu/Run2015B-PromptReco-v1/MINIAOD
#   /SingleMuon/Run2015B-PromptReco-v1/MINIAOD


#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
#config.Data.inputDataset =  '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/MINIAODSIM'
#config.Data.inputDataset =  '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
#config.Data.inputDataset =  '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
#config.Data.inputDataset =  '/WWTo2L2Nu_13TeV-powheg/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'

#config.Data.inputDataset = '/ZZTo4L_Tune4C_13TeV-powheg-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset = '/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset  = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU4bx50_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 4
#config.Data.totalUnits = 25
#config.Data.outLFN = '/store/user/theofil/test' # or '/store/group/<subdir>'
config.Data.outLFNDirBase = '/store/user/theofil/test' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = ''


config.Site.storageSite = 'T2_CH_CERN'

#config.JobType.allowNonProductionCMSSW = True ##  CMSSW_7_3_2 on slc6_amd64_gcc481 is not among supported releases 
