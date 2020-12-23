#ifndef DQM_GEM_INTERFACE_GEMDQMBase_h
#define DQM_GEM_INTERFACE_GEMDQMBase_h

#include <map>
#include <tuple>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "Validation/MuonGEMHits/interface/GEMValidationUtils.h"

#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"


class GEMDQMBase : public DQMEDAnalyzer {
public:
  ~GEMDQMBase() override{};

protected:
  int initGeometry(edm::EventSetup const& iSetup);
  int loadChambers();
  int findVFAT(float min_, float max_, float x_, int roll_);

  int GenerateMEPerChamber(DQMStore::IBooker& ibooker);
  virtual int ProcessWithMEMap2(DQMStore::IBooker& ibooker, ME2IdsKey key) {return 0;};  // must be overrided
  virtual int ProcessWithMEMap3(DQMStore::IBooker& ibooker, ME3IdsKey key) {return 0;};  // must be overrided
  virtual int ProcessWithMEMap4(DQMStore::IBooker& ibooker, ME4IdsKey key) {return 0;};  // must be overrided

  const GEMGeometry* GEMGeometry_;

  std::vector<GEMChamber> gemChambers_;

  std::map<ME2IdsKey, bool> MEMap2Check_;
  std::map<ME3IdsKey, bool> MEMap3Check_;
  std::map<ME4IdsKey, bool> MEMap4Check_;
};

#endif  // DQM_GEM_INTERFACE_GEMDQMBase_h
