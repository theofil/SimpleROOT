import FWCore.ParameterSet.Config as cms

process = cms.Process("demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService", 
      fileName = cms.string("output.root"),
      closeFileFast = cms.untracked.bool(True)
  )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FEE3CF68-796C-E411-ABF5-002590DB9214.root' 
    )
)

process.demo = cms.EDAnalyzer('SimpleROOT')
process.p = cms.Path(process.demo)
