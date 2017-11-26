import FWCore.ParameterSet.Config as cms

GEMDQMSourceDigi = cms.EDAnalyzer("GEMDQMSourceDigi",
    digisInputLabel = cms.InputTag("muonGEMDigis", ""),
    errorsInputLabel = cms.InputTag("muonGEMDigis", "vfatStatus"),
    AMCInputLabel = cms.InputTag("muonGEMDigis", "AMCStatus"),     
    GEBInputLabel = cms.InputTag("muonGEMDigis", "GEBStatus"), 
    
    AMCBID = cms.untracked.vint32(
        -1, -1, 48879, -1, 
        -1, -1, -1, -1, 
        -1, -1, -1, -1, 
    ), 
    
    AMCBinLabel = cms.untracked.vstring(
        "label1", "label2", "48879", "label4", 
        "label5", "label6", "label7", "label8", 
        "label9", "label10", "label11", "label12", 
    )
)
