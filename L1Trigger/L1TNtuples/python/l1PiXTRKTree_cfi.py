import FWCore.ParameterSet.Config as cms

from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import *

l1PiXTRKTree = cms.EDAnalyzer(
    "L1PiXTRKTreeProducer",
    pixelRecHits = cms.InputTag("siPixelRecHits"),
    genParticleToken = cms.InputTag("genParticles"),
    pileupInfoToken     = cms.untracked.InputTag("addPileupInfo"),
    #egCrys = cms.InputTag("L1EGammaCrystalsProducer","EGammaCrystal"),
    egCrys = cms.InputTag("L1EGammaCrystalsProducer","L1EGCollectionWithCuts"),
    egCrysCluster = cms.InputTag("L1EGammaCrystalsProducer","L1EGXtalClusterWithCuts"),
    eg = cms.InputTag("simCaloStage2Digis"), # Run 2

    # HGCal
    TriggerCells = cms.InputTag('hgcalTriggerPrimitiveDigiProducer:calibratedTriggerCells'),
    Clusters = cms.InputTag('hgcalTriggerPrimitiveDigiProducer:cluster2D'),
    Multiclusters = cms.InputTag('hgcalTriggerPrimitiveDigiProducer:cluster3D'),
    
    me0Segment = cms.InputTag("me0TriggerPseudoDigis"), 

    simVertex = cms.InputTag("g4SimHits", "", "SIM"),
    simTrack = cms.InputTag("g4SimHits", "", "SIM"),

    L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
    
    TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
    TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth")

   
)
