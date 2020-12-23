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

#include "DQM/GEM/interface/GEMDQMBase.h"

#include <string>

//----------------------------------------------------------------------------------------------------

//class GEMDigiSource : public DQMEDAnalyzer
class GEMDigiSource : public GEMDQMBase {
public:
  GEMDigiSource(const edm::ParameterSet& cfg);
  ~GEMDigiSource() override{};
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  void dqmBeginRun(edm::Run const&, edm::EventSetup const&) override{};
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(edm::Event const& e, edm::EventSetup const& eSetup) override;

private:
  virtual int ProcessWithMEMap3(DQMStore::IBooker& ibooker, ME3IdsKey key) override;
  
  edm::EDGetToken tagDigi_;

  MEMap3Ids mapStripOcc_ieta_;
  MEMap3Ids mapStripOcc_phi_;
  MEMap3Ids mapTotalStrip_layer_;

  MEMap3Ids mapBX_layer_;
  
  Int_t nBXMin_, nBXMax_;
};

using namespace std;
using namespace edm;

GEMDigiSource::GEMDigiSource(const edm::ParameterSet& cfg) {
  tagDigi_ = consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("digisInputLabel"));
}

void GEMDigiSource::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("digisInputLabel", edm::InputTag("muonGEMDigis", ""));
  descriptions.add("GEMDigiSource", desc);
}

void GEMDigiSource::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const&, edm::EventSetup const& iSetup) {
  initGeometry(iSetup);
  if (GEMGeometry_ == nullptr)
    return;
  loadChambers();

  nBXMin_ = -10;
  nBXMax_ =  10;

  ibooker.cd();
  ibooker.setCurrentFolder("GEM/digi");
  GenerateMEPerChamber(ibooker);
  //for (const auto& ch : gemChambers_) {
  //   FIXME: ge+XY or ge-XY for names, GE+XY or GE-XY for titles
  //  std::string strIdxName = "Gemini_" + to_string(gid.chamber()) + "_GE" + (gid.region() > 0 ? "p" : "m") +
  //                           to_string(gid.station()) + "_" + to_string(gid.layer());
  //  std::string strIdxName = "_GE1" + to_string(gid.station()) +  // "GE1" will be replaced by "GE[12]"
  //         (gid.region() > 0 ? "-P-" : "-M-") + 
  //         to_string(gid.chamber()) + "L" + to_string(gid.layer()) + ( gid.chamber() % 2 == 0 ? "-L" : "-S" );
  //  std::string strIdxTitle = strIdxName;

  //  string hName_digi = "Digi_Strips_" + strIdxName;
  //  string hTitle_digi = "Digi Strip " + strIdxTitle;
  //  string hAxis_digi = ";Strip;iEta";
  //  Digi_2D_[ch.id()] = ibooker.book2D(hName_digi, hTitle_digi + hAxis_digi, 384, 1, 385, 8, 0.5, 8.5);
  //  Digi_1D_[ch.id()] = ibooker.book1D(hName_digi + "_VFAT", hTitle_digi + " VFAT" + hAxis_digi, 24, 0, 24);

  //  string hNameBx = "bx_vs_VFAT_" + strIdxName;
  //  string hTitleBx = "bx vs VFAT " + strIdxTitle;
  //  hTitleBx += ";Bunch crossing;VFAT";
  //  BxVsVFAT[ch.id()] = ibooker.book2D(hNameBx, hTitleBx, 10, -5, 5, 24, 0, 24);
  //}
}

int GEMDigiSource::ProcessWithMEMap3(DQMStore::IBooker& ibooker, ME3IdsKey key) {
  auto strSuffixName  = GEMUtils::getSuffixName(key);
  auto strSuffixTitle = GEMUtils::getSuffixTitle(key);

  mapStripOcc_ieta_[key] = ibooker.book1D(
      "strip_ieta_occ" + strSuffixName, 
      "Strip Digi Occupancy per eta partition" + strSuffixTitle, 
      8, 0.5, 8.5
    );

  mapStripOcc_phi_[key] = ibooker.book1D(
      "strip_phi_occ" + strSuffixName, 
      "Strip Digi Phi (degree) Occupancy " + strSuffixTitle, 
      108, -5, 355
    );

  mapTotalStrip_layer_[key] = ibooker.book1D(
      "total_strips_per_event" + strSuffixName, 
      "Total number of strip digis per event for each layers " + strSuffixTitle, 
      50, -0.5, 99.5
    );

  mapBX_layer_[key] = ibooker.book1D(
      "strip_phi_occ" + strSuffixName, 
      "Strip Digi Bunch Crossing " + strSuffixTitle, 
      11, nBXMin_ - 0.5, nBXMax_ + 0.5
    );

  return 0;
}

void GEMDigiSource::analyze(edm::Event const& event, edm::EventSetup const& eventSetup) {
  edm::Handle<GEMDigiCollection> gemDigis;
  event.getByToken(this->tagDigi_, gemDigis);
  std::map<ME3IdsKey, Int_t> total_strip_layer;
  for (const auto& ch : gemChambers_) {
    GEMDetId gid = ch.id();
    ME2IdsKey key2{gid.region(), gid.station()};
    ME3IdsKey key3{gid.region(), gid.station(), gid.layer()};
    const BoundPlane& surface = GEMGeometry_->idToDet(gid)->surface();
    if (total_strip_layer.find(key3) == total_strip_layer.end()) total_strip_layer[key3] = 0;
    for (auto roll : ch.etaPartitions()) {
      GEMDetId rId = roll->id();
      const auto& digis_in_det = gemDigis->get(rId);
      for (auto d = digis_in_det.first; d != digis_in_det.second; ++d) {
        mapStripOcc_ieta_[key3]->Fill(rId.roll());
        GlobalPoint strip_global_pos = surface.toGlobal(roll->centreOfStrip(d->strip()));
        Float_t fPhiDeg = ( (Float_t)strip_global_pos.phi() ) * 180.0 / 3.141592;
        if (fPhiDeg < -5.0) fPhiDeg += 360.0;
        mapStripOcc_phi_[key3]->Fill(fPhiDeg);
        total_strip_layer[key3]++;
        Int_t nBX = std::min(std::max((Int_t)d->bx(), nBXMin_), nBXMax_);
        mapBX_layer_[key3]->Fill(nBX);  // Over/underflow
      }
    }
  }
  for (auto [key, num_total_strip] : total_strip_layer) mapTotalStrip_layer_[key]->Fill(num_total_strip);
}

DEFINE_FWK_MODULE(GEMDigiSource);
