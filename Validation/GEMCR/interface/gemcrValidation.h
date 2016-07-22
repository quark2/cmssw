#ifndef gemcrValidation_H
#define gemcrValidation_H

#include "Validation/MuonGEMHits/interface/GEMBaseValidation.h"
#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include <DataFormats/GEMRecHit/interface/GEMRecHitCollection.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMChamber.h"


class gemcrValidation : public GEMBaseValidation
{
public:
  explicit gemcrValidation( const edm::ParameterSet& );
  ~gemcrValidation();
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event& e, const edm::EventSetup&) override;
  int findIndex(GEMDetId id_);
  const GEMGeometry* initGeometry(edm::EventSetup const & iSetup);
  //MonitorElement* BookHist1D( DQMStore::IBooker &, const char* name, const char* label, unsigned int region_num, unsigned int station_num, unsigned int layer_num, const unsigned int Nbin, const Float_t xMin, const Float_t xMax);
  //MonitorElement* BookHist1D( DQMStore::IBooker &, const char* name, const char* label, unsigned int region_num, const unsigned int Nbin, const Float_t xMin, const Float_t xMax);

private:
  const GEMGeometry* GEMGeometry_;
  //Detaile Plots
  //MonitorElement* gem_rh_xy[2][3][2];
  //MonitorElement* gem_rh_zr[2][3][2];
  //MonitorElement* gem_cls[2][3][2];
  //MonitorElement* gem_pullX[2][3][2];
  //MonitorElement* gem_pullY[2][3][2];
  //std::vector<MonitorElement*> gem_chamber_x_y;
  std::vector<MonitorElement*> gem_chamber_x_y;
  std::vector<MonitorElement*> gem_chamber_cl_size;
  std::vector<MonitorElement*> gem_chamber_bx;
  std::vector<MonitorElement*> gem_chamber_firedStrip;
  std::vector<MonitorElement*> gem_chamber_tr_eff;
  std::vector<MonitorElement*> gem_chamber_ch_eff;
  
  MonitorElement* gemcr_g;
  
  //Simple Plots
  MonitorElement* gem_cls_tot;
  MonitorElement* gem_bx_tot;
  MonitorElement* tr_eff;
  MonitorElement* hit_eff;
  MonitorElement* tr_size;
  MonitorElement* tr_hit_size;
  MonitorElement* tr_chamber;

  std::vector<GEMChamber> gemChambers;
  int n_ch;
  MuonServiceProxy* theService;
 

  //std::unordered_map< UInt_t , MonitorElement* > recHits_dcEta;
  //std::unordered_map< UInt_t , MonitorElement* > recHits_simple_zr;
  //MonitorElement* gem_region_pullX[2];
  //MonitorElement* gem_region_pullY[2];

  edm::EDGetToken InputTagToken_, InputTagToken_RH, InputTagToken_TR, InputTagToken_TS;

};

#endif
