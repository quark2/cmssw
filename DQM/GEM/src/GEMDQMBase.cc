#include "DQM/GEM/interface/GEMDQMBase.h"


using namespace std;
using namespace edm;


//GEMDQMBase::GEMDQMBase(const edm::ParameterSet& cfg) {}


int GEMDQMBase::initGeometry(edm::EventSetup const& iSetup) {
  GEMGeometry_ = nullptr;
  try {
    edm::ESHandle<GEMGeometry> hGeom;
    iSetup.get<MuonGeometryRecord>().get(hGeom);
    GEMGeometry_ = &*hGeom;
  } catch (edm::eventsetup::NoProxyException<GEMGeometry>& e) {
    edm::LogError("MuonGEMBaseValidation") << "+++ Error : GEM geometry is unavailable on event loop. +++\n";
    return -1;
  }

  return 0;
}


int GEMDQMBase::loadChambers() {
  if ( GEMGeometry_ == nullptr ) return -1;
  gemChambers_.clear();
  const std::vector<const GEMSuperChamber*>& superChambers_ = GEMGeometry_->superChambers();
  for (auto sch : superChambers_) {
    int n_lay = sch->nChambers();
    for (int l = 0; l < n_lay; l++) {
      Bool_t bExist = false;
      for (const auto& ch : gemChambers_) {
        if (ch.id() == sch->chamber(l + 1)->id()) {
          bExist = true;
          break;
        }
      }
      if (bExist) continue;
      gemChambers_.push_back(*sch->chamber(l + 1));
    }
  }
  return 0;
}


int GEMDQMBase::findVFAT(float min_, float max_, float x_, int roll_) {
  float step = max_ / 3;
  if (x_ < (min_ + step)) {
    return 8 - roll_;
  } else if (x_ < (min_ + 2 * step)) {
    return 16 - roll_;
  } else {
    return 24 - roll_;
  }
}


int GEMDQMBase::GenerateMEPerChamber(DQMStore::IBooker& ibooker) {
  MEMap2Check_.clear();
  MEMap3Check_.clear();
  MEMap4Check_.clear();
  for (const auto& ch : gemChambers_) {
    GEMDetId gid = ch.id();
    ME2IdsKey key2{gid.region(), gid.station()};
    ME3IdsKey key3{gid.region(), gid.station(), gid.layer()};
    if (!MEMap2Check_[ key2 ]) {
      ProcessWithMEMap2(ibooker, key2);
      MEMap2Check_[ key2 ] = true;
    }
    if (!MEMap3Check_[ key3 ]) {
      ProcessWithMEMap3(ibooker, key3);
      MEMap3Check_[ key3 ] = true;
    }
    for (auto roll : ch.etaPartitions()) {
      GEMDetId rId = roll->id();
      ME4IdsKey key4{gid.region(), gid.station(), gid.layer(), rId.roll()};
      if (!MEMap4Check_[ key4 ]) {
        ProcessWithMEMap4(ibooker, key4);
        MEMap4Check_[ key4 ] = true;
      }
    }
  }
  return 0;
}


