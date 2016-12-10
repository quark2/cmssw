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

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include <TFile.h>
#include <TEfficiency.h>
#include <TH1D.h>

///
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/GEMGeometry/interface/GEMGeometry.h>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/CosmicMuonProducer/interface/CosmicMuonSmoother.h"
#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneMuonSmoother.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
///









class gemcrValidation : public GEMBaseValidation
{
public:
  explicit gemcrValidation( const edm::ParameterSet& );
  ~gemcrValidation();
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event& e, const edm::EventSetup&) override;
  int findIndex(GEMDetId id_);
  int findvfat(float x, float a, float b);
  const GEMGeometry* initGeometry(edm::EventSetup const & iSetup);
  double maxCLS, minCLS,maxRes;
private:
  const GEMGeometry* GEMGeometry_;
  std::vector<MonitorElement*> gem_chamber_x_y;
  std::vector<MonitorElement*> gem_chamber_cl_size;
  std::vector<MonitorElement*> gem_chamber_bx;
  std::vector<MonitorElement*> gem_chamber_firedStrip;
  std::vector<MonitorElement*> gem_chamber_tr2D_eff;
  std::vector<MonitorElement*> gem_chamber_th2D_eff;
  std::vector<MonitorElement*> gem_chamber_trxy_eff;
  std::vector<MonitorElement*> gem_chamber_thxy_eff;
  std::vector<MonitorElement*> gem_chamber_residual;
  std::vector<MonitorElement*> gem_chamber_residual_r;
  std::vector<MonitorElement*> gem_chamber_local_x;
  
  MonitorElement* gemcr_g;
  MonitorElement* gem_cls_tot;
  MonitorElement* gem_bx_tot;
  MonitorElement* tr_size;
  MonitorElement* tr_hit_size;
  MonitorElement* tr_chamber;
  MonitorElement* th_chamber;
  MonitorElement* rh_chamber;
  MonitorElement* rh1_chamber;
  MonitorElement* rh2_chamber;
  MonitorElement* rh3_chamber;

  std::vector<GEMChamber> gemChambers;
  int n_ch;
  
  MuonServiceProxy* theService;
  
  CosmicMuonSmoother* theSmoother;
  KFUpdator* theUpdator;
  std::auto_ptr<std::vector<TrajectorySeed> > findSeeds(MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits);
  Trajectory makeTrajectory(TrajectorySeed seed, MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits, std::vector<GEMChamber> gemChambers, GEMChamber testChamber);
  edm::EDGetToken InputTagToken_, InputTagToken_RH, InputTagToken_TR, InputTagToken_TS, InputTagToken_GP;
};

#endif
