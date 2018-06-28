import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

from L1Trigger.L1TNtuples.l1PiXTRKTree_cfi import *

l1PiXTRKTree.ecalToken = cms.untracked.InputTag("simEcalTriggerPrimitiveDigis")
l1PiXTRKTree.hcalToken = cms.untracked.InputTag("simHcalTriggerPrimitiveDigis")
l1PiXTRKTree.l1TowerToken = cms.untracked.InputTag("simCaloStage2Layer1Digis")
l1PiXTRKTree.l1ClusterToken = cms.untracked.InputTag("simCaloStage2Digis", "MP")

L1NtuplePiXTRK = cms.Sequence(
  l1PiXTRKTree
)
