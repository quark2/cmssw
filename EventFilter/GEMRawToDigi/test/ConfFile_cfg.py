import FWCore.ParameterSet.Config as cms

process = cms.Process("RAW")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1200)) # 370000

process.source = cms.Source("EmptySource"
#    setRunNumber=cms.untracked.uint32(999)
)


process.gemDigis = cms.EDProducer('GEMFedRawProducerV2',
      inputFileName=cms.string("GEMDQMRawData.dat"),
      slotVector=cms.vint32(0x0,0xeac,0x0,0x0,0x0,0x0,0xe67,0xea0,0x0,0x0,0x0,0x0,0x0,0xe87,0x0,0x0,0x0,0x0,0xe8f,0xe93,0x0,0x0,0xe84,0xe88),
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myGEMFedRawDataFile.root')
)

  
process.p = cms.Path(process.gemDigis)

process.e = cms.EndPath(process.out)
