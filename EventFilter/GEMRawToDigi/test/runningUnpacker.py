import FWCore.ParameterSet.Config as cms

process = cms.Process("RAW")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000))

process.source = cms.Source("EmptySource"   # This should probably be changed to another kind of source
#    setRunNumber=cms.untracked.uint32(999)
)

process.load("EventFilter.GEMRawToDigi.GEMSQLiteCabling_cfi")
process.GEMCabling.connect = 'sqlite_file:GEMEMap_CosmicStand_7Nov2016.db'
#process.GEMCabling.connect = 'sqlite_file:GEMEMap.db'

SLOT=  [20,    21,    22]
VFAT=  [0xfa28,0xfa17,0xfa17]
COLUMN=[1,     1,     1]
ROW=   [3,     3,     3]
LAYER= [1,     1,     1]

process.gemDigis = cms.EDProducer('GEMCosmicStandUnpacker',
    inputFileName=cms.string("/afs/cern.ch/work/c/cepeda/performance/GEM/CMSSW_8_1_0_pre12/src/datatestJared/run000047_teststand_CERN904_2016-09-26_chunk_99.dat"),
#      inputFileName=cms.string("GEMDQMRawData.dat"),
      slotVector=  cms.vint32(SLOT),
      vfatVector=cms.vuint64(VFAT),
      columnVector=cms.vint32(COLUMN),
      rowVector   =cms.vint32(ROW),
      layerVector = cms.vint32(LAYER)
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myGEMFedRawDataFile.root')
)

process.check = cms.EDAnalyzer("GEMDigiReader",
      InputLabel = cms.InputTag("gemDigis"))


  
process.p = cms.Path(process.gemDigis)#*process.check)

process.e = cms.EndPath(process.out)
