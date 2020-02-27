import FWCore.ParameterSet.Config as cms

import os

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('DQMTEST')


process.MessageLogger = cms.Service("MessageLogger",
  statistics = cms.untracked.vstring(),
  destinations = cms.untracked.vstring('cerr'),
  cerr = cms.untracked.PSet(
      threshold = cms.untracked.string('WARNING')
  )
)

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
#process.load('Configuration.Geometry.GeometryExtended2026D35Reco_cff')
process.load("DQM.Integration.config.FrontierCondition_GT_cfi")


process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    'file:/eos/cms/store/data/Commissioning2020/Cosmics/RAW/v1/000/335/685/00000/2DC31764-CAFB-9946-91FD-56452139E982.root'
  ),
  inputCommands = cms.untracked.vstring(
    'keep *',
  )
)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

process.load("EventFilter.GEMRawToDigi.muonGEMDigis_cfi")
process.load('RecoLocalMuon.GEMRecHit.gemRecHits_cfi')

# dump raw data
process.dumpRaw = cms.EDAnalyzer(
    "DumpFEDRawDataProduct",
    token = cms.untracked.InputTag("rawDataCollector"),
    feds = cms.untracked.vint32 ( 1467,1468 ),
    dumpPayload = cms.untracked.bool ( False )
)

process.output = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring("keep *"),
    #fileName = cms.untracked.string('gem_EDM.root')
    fileName = cms.untracked.string('gem_EDM.root')
)

############## DB file ################# 
#from CondCore.CondDB.CondDB_cfi import *
#CondDB.DBParameters.authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
#CondDB.connect = cms.string('sqlite_fip:DQM/GEM/data/GEMeMap.db')
#
#process.GEMCabling = cms.ESSource("PoolDBESSource",
#    CondDB,
#    toGet = cms.VPSet(cms.PSet(
#        record = cms.string('GEMeMapRcd'),
#        tag = cms.string('GEMeMap_v2')
#    )),
#)
####################################
process.path = cms.Path(
  process.muonGEMDigis *
  process.gemRecHits
)

process.end_path = cms.EndPath(
  process.output
)

process.schedule = cms.Schedule(
  process.path,
  process.end_path
)
