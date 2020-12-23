//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/PluginManager/interface/ModuleDef.h"
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Utilities/interface/InputTag.h"
//#include "FWCore/MessageLogger/interface/MessageLogger.h"
//
//#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
//#include "DQMServices/Core/interface/DQMStore.h"
//#include "DQMServices/Core/interface/MonitorElement.h"
//
#include "DataFormats/GEMRecHit/interface/GEMRecHit.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
//
//#include "Validation/MuonGEMHits/interface/GEMValidationUtils.h"
//
//#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
//#include "Geometry/Records/interface/MuonGeometryRecord.h"
//
//#include <string>
#include "DQM/GEM/interface/GEMDQMBase.h"

#include <string>

//----------------------------------------------------------------------------------------------------

class GEMRecHitSource : public GEMDQMBase {
public:
  GEMRecHitSource(const edm::ParameterSet& cfg);
  ~GEMRecHitSource() override{};
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  void dqmBeginRun(edm::Run const&, edm::EventSetup const&) override{};
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(edm::Event const& e, edm::EventSetup const& eSetup) override;

private:
  virtual int ProcessWithMEMap3(DQMStore::IBooker& ibooker, ME3IdsKey key) override;
  virtual int ProcessWithMEMap4(DQMStore::IBooker& ibooker, ME4IdsKey key) override;

  edm::EDGetToken tagRecHit_;

  float fGlobXMin_, fGlobXMax_;
  float fGlobYMin_, fGlobYMax_;

  int nIdxFirstStrip_;
  int nClusterSizeBinNum_;

  MEMap3Ids mapRecHitOcc_ieta_;
  MEMap3Ids mapRecHitOcc_phi_;
  MEMap3Ids mapTotalRecHit_layer_;
  MEMap4Ids mapCLSRecHit_roll_;

  Int_t nCLSMax_;

  std::unordered_map<UInt_t, MonitorElement*> recHitME_;
  std::unordered_map<UInt_t, MonitorElement*> VFAT_vs_ClusterSize_;
  std::unordered_map<UInt_t, MonitorElement*> StripsFired_vs_eta_;
  std::unordered_map<UInt_t, MonitorElement*> rh_vs_eta_;
  std::unordered_map<UInt_t, MonitorElement*> recGlobalPos;
};

using namespace std;
using namespace edm;

GEMRecHitSource::GEMRecHitSource(const edm::ParameterSet& cfg) {
  tagRecHit_ = consumes<GEMRecHitCollection>(cfg.getParameter<edm::InputTag>("recHitsInputLabel"));

  nIdxFirstStrip_ = cfg.getParameter<int>("idxFirstStrip");
  nClusterSizeBinNum_ = cfg.getParameter<int>("ClusterSizeBinNum");

  fGlobXMin_ = cfg.getParameter<double>("global_x_bound_min");
  fGlobXMax_ = cfg.getParameter<double>("global_x_bound_max");
  fGlobYMin_ = cfg.getParameter<double>("global_y_bound_min");
  fGlobYMax_ = cfg.getParameter<double>("global_y_bound_max");
}

void GEMRecHitSource::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("recHitsInputLabel", edm::InputTag("gemRecHits", ""));

  desc.add<int>("idxFirstStrip", 0);
  desc.add<int>("ClusterSizeBinNum", 9);

  desc.add<double>("global_x_bound_min", -350);
  desc.add<double>("global_x_bound_max", 350);
  desc.add<double>("global_y_bound_min", -260);
  desc.add<double>("global_y_bound_max", 260);

  descriptions.add("GEMRecHitSource", desc);
}

void GEMRecHitSource::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const&, edm::EventSetup const& iSetup) {
  std::vector<GEMDetId> listLayerOcc;

  initGeometry(iSetup);
  if (GEMGeometry_ == nullptr)
    return;
  loadChambers();

  ibooker.cd();
  ibooker.setCurrentFolder("GEM/recHit");
  
  nCLSMax_ = 10;

  GenerateMEPerChamber(ibooker);
}

int GEMRecHitSource::ProcessWithMEMap3(DQMStore::IBooker& ibooker, ME3IdsKey key) {
  auto strSuffixName  = GEMUtils::getSuffixName(key);
  auto strSuffixTitle = GEMUtils::getSuffixTitle(key);

  mapRecHitOcc_ieta_[key] = ibooker.book1D(
      "rechit_ieta_occ" + strSuffixName, 
      "RecHit Digi Occupancy per eta partition" + strSuffixTitle, 
      8, 0.5, 8.5
    );

  mapRecHitOcc_phi_[key] = ibooker.book1D(
      "rechit_phi_occ" + strSuffixName, 
      "RecHit Digi Phi (degree) Occupancy " + strSuffixTitle, 
      108, -5, 355
    );

  mapTotalRecHit_layer_[key] = ibooker.book1D(
      "total_rechit_per_event" + strSuffixName, 
      "Total number of rec. hits per event for each layers " + strSuffixTitle, 
      50, -0.5, 99.5
    );

  return 0;
}

int GEMRecHitSource::ProcessWithMEMap4(DQMStore::IBooker& ibooker, ME4IdsKey key) {
  auto strSuffixName  = GEMUtils::getSuffixName(key);
  auto strSuffixTitle = GEMUtils::getSuffixTitle(key);

  mapCLSRecHit_roll_[key] = ibooker.book1D(
      "total_rechit_per_event" + strSuffixName, 
      "Total number of rec. hits per event for each layers " + strSuffixTitle, 
      nCLSMax_, -0.5, nCLSMax_ + 0.5
    );

  return 0;
}

void GEMRecHitSource::analyze(edm::Event const& event, edm::EventSetup const& eventSetup) {
  edm::Handle<GEMRecHitCollection> gemRecHits;
  event.getByToken(this->tagRecHit_, gemRecHits);
  if (!gemRecHits.isValid()) {
    edm::LogError("GEMRecHitSource") << "GEM recHit is not valid.\n";
    return;
  }
  std::map<ME3IdsKey, Int_t> total_rechit_layer;
  for (const auto& ch : gemChambers_) {
    GEMDetId gid = ch.id();
    ME3IdsKey key3{gid.region(), gid.station(), gid.layer()};
    for (auto roll : ch.etaPartitions()) {
      GEMDetId rId = roll->id();
      ME4IdsKey key4{gid.region(), gid.station(), gid.layer(), rId.roll()};
      if (total_rechit_layer.find(key3) == total_rechit_layer.end()) total_rechit_layer[key3] = 0;
      const auto& recHitsRange = gemRecHits->get(rId);
      auto gemRecHit = recHitsRange.first;
      for (auto hit = gemRecHit; hit != recHitsRange.second; ++hit) {
        mapRecHitOcc_ieta_[key3]->Fill(rId.roll());
        GlobalPoint recHitGP = GEMGeometry_->idToDet(hit->gemId())->surface().toGlobal(hit->localPosition());
        Float_t fPhiDeg = ( (Float_t)recHitGP.phi() ) * 180.0 / 3.141592;
        if (fPhiDeg < -5.0) fPhiDeg += 360.0;
        mapRecHitOcc_phi_[key3]->Fill(fPhiDeg);
        total_rechit_layer[key3]++;
        mapCLSRecHit_roll_[key4]->Fill(std::min(hit->clusterSize(), nCLSMax_));
      }
    }
  }
  for (auto [key, num_total_rechit] : total_rechit_layer) mapTotalRecHit_layer_[key]->Fill(num_total_rechit);
}

DEFINE_FWK_MODULE(GEMRecHitSource);
