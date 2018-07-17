process = cms.Process("PickEvent")

process.source = cms.Source ("PoolSource",
  #fileNames = cms.untracked.vstring("file:/xrootd/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/4E043DF2-A4BE-E611-84AF-0025905B85C0.root")
  #fileNames = cms.untracked.vstring("file:4E043DF2-A4BE-E611-84AF-0025905B85C0.root")
  fileNames = cms.untracked.vstring("file:/xrootd/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/C2EBB5C9-1AED-E611-814B-0CC47A7E6B00.root")
)

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32 (10000))

process.metAnalyzer = cms.EDAnalyzer('METAnalyzer',
  src = cms.InputTag('slimmedMETs'), 
)

process.run = cms.Path(process.metAnalyzer)

#process.Out = cms.OutputModule("PoolOutputModule",
#  fileName = cms.untracked.string ("MyOutputFile.root")
#)
process.TFileService = cms.Service("TFileService", fileName = cms.string("out.root"))

#process.end = cms.EndPath(process.Out)
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
