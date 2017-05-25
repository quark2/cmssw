import FWCore.ParameterSet.Config as cms

from DQM.GEM.GEMDQMSource_cfi import *

GEMDQM = cms.Sequence(
  GEMDQMSource
)
