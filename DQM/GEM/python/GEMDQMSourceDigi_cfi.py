import FWCore.ParameterSet.Config as cms

GEMDQMSourceDigi = cms.EDAnalyzer("GEMDQMSourceDigi",
    digisInputLabel = cms.InputTag("muonGEMDigis", ""),
    errorsInputLabel = cms.InputTag("muonGEMDigis", "vfatStatus"),
    AMCInputLabel = cms.InputTag("muonGEMDigis", "AMCStatus"),     
    GEBInputLabel = cms.InputTag("muonGEMDigis", "GEBStatus"), 
)
