import FWCore.ParameterSet.Config as cms

muonGEMDigis = cms.EDProducer("GEMRawToDigiModule",
    InputLabel = cms.InputTag("rawDataCollector"),
    UnpackStatusDigis = cms.untracked.bool(True)
    useDBEMap = cms.bool(True),    
)
