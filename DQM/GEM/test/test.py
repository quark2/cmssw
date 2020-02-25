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
#process.GlobalTag.globaltag = '110X_dataRun3_HLT_v1'



process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "GEM"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "GEM"

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    #'file:/eos/cms/store/express/Commissioning2018/ExpressCosmics/FEVT/Express-v1/000/310/292/00000/6C23251D-4F18-E811-AEC5-02163E01A41D.root'
    #'file:/cms/ldap_home/quark2930/Work/dqm/dqmvalidation_latest_200123/src/20011.0_new/step3.root'
    'file:/eos/cms/store/data/Commissioning2020/Cosmics/RAW/v1/000/335/685/00000/2DC31764-CAFB-9946-91FD-56452139E982.root'
  ),
  inputCommands = cms.untracked.vstring(
    'keep *',
    #'keep FEDRawDataCollection_*_*_*'
  )
)

import sys

strNameList = sys.argv[ 2 ] if len(sys.argv) > 2 else ""

if strNameList != "": 
  listSrc = []
  if strNameList.endswith(".root"): 
    listSrc.append(strNameList)
  else: 
    with open(strNameList) as fSrc: listSrc += [ s for s in fSrc.read().splitlines() if s != "" ]
  
  process.source.fileNames = cms.untracked.vstring(listSrc)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

process.load("EventFilter.GEMRawToDigi.muonGEMDigis_cfi")
process.load('RecoLocalMuon.GEMRecHit.gemRecHits_cfi')
process.load("DQM.GEM.GEMDQM_cff")

process.muonGEMDigis.useDBEMap = True
process.muonGEMDigis.unPackStatusDigis = True

#process.GEMDQMStatusDigi.pathOfPrevDQMRoot = "DQM_V0001_GEM_R000030000.root"
#process.GEMDQMHarvester.fromFile = "DQM_V0001_GEM_R000020150.root"

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
  process.gemRecHits *
  process.GEMDQM
)

process.end_path = cms.EndPath(
  process.dqmEnv +
  process.dqmSaver
)

process.schedule = cms.Schedule(
  process.path,
  process.end_path
)
