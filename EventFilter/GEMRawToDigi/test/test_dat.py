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



process.source = cms.Source(
    "NewEventStreamFileReader",
    fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/b/bko/Public/dat_run335685/run335685_ls0001_streamDQM_mrg-c2f12-23-01.dat"),
    skipEvents=cms.untracked.uint32(0)
)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

process.load("EventFilter.GEMRawToDigi.muonGEMDigis_cfi")
process.load('RecoLocalMuon.GEMRecHit.gemRecHits_cfi')

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
