import FWCore.ParameterSet.Config as cms

process = cms.Process("RAW")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000))

process.source = cms.Source("EmptySource",   # This should probably be changed to another kind of source
      firstEvent = cms.untracked.uint32(1),
      firstRun = cms.untracked.uint32(1)
)

process.load("EventFilter.GEMRawToDigi.GEMSQLiteCabling_cfi")
process.GEMCabling.connect = 'sqlite_file:GEMEMap_CosmicStand_7Nov2016.db'
#process.GEMCabling.connect = 'sqlite_file:GEMEMap.db'

VFAT=[0xf943, 0xf743, 0xf787, 0xf6cb, 0xffdf, 0xff44, 0xf6c4, 0xfa13, 0xffc4, 0xf6b0, 0xf6cc, 0xfe27, 0xf920, 0xf968, 0xff90, 0xff3c, 0xf783, 0xff67, 0xffc7, 0xf6dc, 0xff20, 0xf6a4, 0xff58, 0xffd7, 0xf94b, 0xf6b4, 0xf6a3, 0xf72f, 0xfa07, 0xff9c, 0xfec8, 0xffdb, 0xfe33, 0xf6e3, 0xff64, 0xfff0, 0xfe43, 0xf98f, 0xf784, 0xf940, 0xff54, 0xff17, 0xf6e8, 0xf78b, 0xf9e8, 0xf6ff, 0xf958, 0xff50]
SLOT=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,23]
COLUMN=[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
ROW=[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
LAYER=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

process.gemDigis = cms.EDProducer('GEMCosmicStandUnpacker',
      inputFileName=cms.string("/data/bigdisk/GEM-Data-Taking/GE11_QC8/run000001_Cosmics_TIF_2016-11-17_chunk_457.dat"),
#    inputFileName=cms.string("/afs/cern.ch/work/c/cepeda/performance/GEM/CMSSW_8_1_0_pre12/src/datatestJared/run000047_teststand_CERN904_2016-09-26_chunk_99.dat"),
#      inputFileName=cms.string("GEMDQMRawData.dat"),
      slotVector=  cms.vint32(SLOT),
      vfatVector=cms.vuint64(VFAT),
      columnVector=cms.vint32(COLUMN),
      rowVector   =cms.vint32(ROW),
      layerVector = cms.vint32(LAYER),
#      verbose=cms.untracked.bool(True)
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myGEMFedRawDataFile.root')
)

process.check = cms.EDAnalyzer("GEMDigiReader",
      InputLabel = cms.InputTag("gemDigis"))


  
process.p = cms.Path(process.gemDigis)#*process.check)

process.e = cms.EndPath(process.out)
