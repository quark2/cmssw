import FWCore.ParameterSet.Config as cms

GEMDQMStatusDigi = cms.EDAnalyzer("GEMDQMStatusDigi",
    VFATInputLabel = cms.InputTag("muonGEMDigis", "vfatStatus"),
    AMCInputLabel = cms.InputTag("muonGEMDigis", "AMCStatus"),     
    GEBInputLabel = cms.InputTag("muonGEMDigis", "GEBStatus"), 
)
