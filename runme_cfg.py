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
   #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/00759690-D16E-E511-B29E-00261894382D.root'
   #'/store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v4/000/258/159/00000/6CA1C627-246C-E511-8A6A-02163E014147.root'
   '/store/data/Run2015D/MET/MINIAOD/PromptReco-v4/000/258/159/00000/1E5A2F7F-D16B-E511-9AC0-02163E0135AC.root'
    )
)

process.demo = cms.EDAnalyzer('SimpleROOT')
process.p = cms.Path(process.demo)
