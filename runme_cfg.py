import FWCore.ParameterSet.Config as cms

process = cms.Process("demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService", 
      fileName = cms.string("output.root"),
      closeFileFast = cms.untracked.bool(True)
  )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/EC2F8972-6C08-E511-B18C-0025905A6090.root'
#       '/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/638/00000/327469FD-0D2B-E511-96BD-02163E012B30.root'
#       '/store/relval/CMSSW_7_4_1/RelValTTbar_13/MINIAODSIM/PU25ns_MCRUN2_74_V9_gensim71X-v1/00000/72C84BA7-F9EC-E411-B875-002618FDA210.root'
#       '/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/FC4E6E16-5C7F-E411-8843-002590200AE4.root',
       #'/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/08B36E8F-5E7F-E411-9D5A-002590200AE4.root'
       #'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FEE3CF68-796C-E411-ABF5-002590DB9214.root'
# '/store/mc/Phys14DR/ZZTo4L_Tune4C_13TeV-powheg-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FA1604CC-E269-E411-9D7F-848F69FD2970.root',
# '/store/mc/Phys14DR/ZZTo4L_Tune4C_13TeV-powheg-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F07DA573-E269-E411-A76A-848F69FD290A.root'
    )
)

process.demo = cms.EDAnalyzer('SimpleROOT')
process.p = cms.Path(process.demo)
