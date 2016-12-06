#include "Validation/GEMCR/interface/gemcrValidation.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "Geometry/GEMGeometry/interface/GEMSuperChamber.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GeometrySurface/interface/Bounds.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Common/interface/RefToBase.h" 
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//#include "RecoMuon/CosmicMuonProducer/interface/CosmicMuonSmoother.h"

//////
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
//#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneMuonSmoother.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
////




#include <iomanip>

#include <TCanvas.h>
using namespace std;

gemcrValidation::gemcrValidation(const edm::ParameterSet& cfg): GEMBaseValidation(cfg)
{
  InputTagToken_RH = consumes<GEMRecHitCollection>(cfg.getParameter<edm::InputTag>("recHitsInputLabel"));
  InputTagToken_TR = consumes<vector<reco::Track>>(cfg.getParameter<edm::InputTag>("tracksInputLabel"));
  InputTagToken_TS = consumes<vector<TrajectorySeed>>(cfg.getParameter<edm::InputTag>("seedInputLabel"));
  edm::ParameterSet serviceParameters = cfg.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  minCLS = cfg.getParameter<double>("minClusterSize"); 
  maxCLS = cfg.getParameter<double>("maxClusterSize");
  maxRes = cfg.getParameter<double>("maxResidual");
  //
  edm::ParameterSet smootherPSet = cfg.getParameter<edm::ParameterSet>("MuonSmootherParameters");
  theSmoother = new CosmicMuonSmoother(smootherPSet,theService);
  theUpdator = new KFUpdator();


}

void gemcrValidation::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & Run, edm::EventSetup const & iSetup ) {
  GEMGeometry_ = initGeometry(iSetup);
  if ( GEMGeometry_ == nullptr) return ;  

  const std::vector<const GEMSuperChamber*>& superChambers_ = GEMGeometry_->superChambers();   
  for (auto sch : superChambers_){
    int n_lay = sch->nChambers();
    for (int l=0;l<n_lay;l++){
      gemChambers.push_back(*sch->chamber(l+1));
    }
  }
  n_ch = gemChambers.size();
  
  ibooker.setCurrentFolder("MuonGEMRecHitsV/GEMRecHitsTask");

  gemcr_g = ibooker.book3D("gemcr_g","GEMCR GLOBAL RECHITS", 200,-100,100,20,-55,55, 140,0,140);
  gem_cls_tot = ibooker.book1D("cluster_size","Cluseter size",20,0,20);  
  gem_bx_tot = ibooker.book1D("bx", "BX" , 30, -15,15);
  tr_size = ibooker.book1D("tr_size", "track size",10,0,10);
  tr_hit_size = ibooker.book1D("tr_hit_size", "hit size in track",15,0,15); 

  tr_chamber = ibooker.book1D("tr_eff_ch", "tr rec /chamber",n_ch,0,n_ch); 
  th_chamber = ibooker.book1D("th_eff_ch", "tr hit/chamber",n_ch,0,n_ch); 
  rh_chamber = ibooker.book1D("rh_eff_ch", "rec hit/chamber",n_ch,0,n_ch); 
  rh1_chamber = ibooker.book1D("rh1_chamber", "all recHits",n_ch,0,n_ch); 
  rh2_chamber = ibooker.book1D("rh2_chamber", "cut passed recHits ",n_ch,0,n_ch); 
  rh3_chamber = ibooker.book1D("rh3_chamber", "tracking recHits",n_ch,0,n_ch); 
  for(int c = 0; c<n_ch; c++){
   cout << gemChambers[c].id() << endl;
   GEMDetId gid = gemChambers[c].id();
   string b_name = "chamber_"+to_string(gid.chamber())+"_layer_"+to_string(gid.layer());
   tr_chamber->setBinLabel(c+1,b_name);
   th_chamber->setBinLabel(c+1,b_name);
   rh_chamber->setBinLabel(c+1,b_name);
   rh1_chamber->setBinLabel(c+1,b_name);
   rh2_chamber->setBinLabel(c+1,b_name);
   rh3_chamber->setBinLabel(c+1,b_name);
  }
  for(int c = 0; c<n_ch;c++){
     GEMDetId gid = gemChambers[c].id();
     string h_name = "chamber_"+to_string(gid.chamber())+"_layer_"+to_string(gid.layer());
     gem_chamber_x_y.push_back(ibooker.book2D(h_name+"_recHit",h_name+" recHit", 500,-25,25,8,1,9));
     gem_chamber_cl_size.push_back(ibooker.book2D(h_name+"_recHit_size", h_name+" recHit size", 50,0,50,24,0,24));
     gem_chamber_firedStrip.push_back(ibooker.book2D(h_name+"_gemDigi", h_name+" gemDigi", 385,0,385,8,1,9));
     gem_chamber_bx.push_back(ibooker.book2D(h_name+"_bx", h_name+" BX", 30,-15,15,10,0,10));
     gem_chamber_tr2D_eff.push_back(ibooker.book2D(h_name+"_recHit_efficiency", h_name+" recHit efficiency", 3,1,4,8,1,9));
     gem_chamber_th2D_eff.push_back(ibooker.book2D(h_name+"_th2D_eff", h_name+"_th2D_eff", 3,1,4,8,1,9));
     gem_chamber_trxy_eff.push_back(ibooker.book2D(h_name+"_trxy_eff", h_name+" recHit efficiency", 50,-25,25,8,1,9));
     gem_chamber_thxy_eff.push_back(ibooker.book2D(h_name+"_thxy_eff", h_name+"_th2D_eff", 50,-25,25,8,1,9));
     //gem_chamber_residual.push_back(ibooker.book2D(h_name+"_residual", h_name+" residual", 140,-7,7,400,-20,20));
     gem_chamber_residual.push_back(ibooker.book2D(h_name+"_residual", h_name+" residual", 500,-25,25,100,-5,5));
     gem_chamber_residual_r.push_back(ibooker.book1D(h_name+"_residual_r", h_name+" residual R", 140,-7,7));
     gem_chamber_local_x.push_back(ibooker.book2D(h_name+"_local_x", "LocalX vs Det_N_LocalX",500,-25,25,500,-25,25));
  }
}

int gemcrValidation::findIndex(GEMDetId id_) {
  int index=-1;
  for(int c =0;c<n_ch;c++){
    if((gemChambers[c].id().chamber() == id_.chamber())&(gemChambers[c].id().layer() == id_.layer()) ){index = c;}
  }
  return index;
}

int gemcrValidation::findvfat(float x, float a, float b) {
  float step = abs(b-a)/3.0;
  if ( x < (min(a,b)+step) ) { return 1;}
  else if ( x < (min(a,b)+2.0*step) ) { return 2;}
  else { return 3;}
}

const GEMGeometry* gemcrValidation::initGeometry(edm::EventSetup const & iSetup) {
  const GEMGeometry* GEMGeometry_ = nullptr;
  try {
    edm::ESHandle<GEMGeometry> hGeom;
    iSetup.get<MuonGeometryRecord>().get(hGeom);
    GEMGeometry_ = &*hGeom;
  }
  catch( edm::eventsetup::NoProxyException<GEMGeometry>& e) {
    edm::LogError("MuonGEMBaseValidation") << "+++ Error : GEM geometry is unavailable on event loop. +++\n";
    return nullptr;
  }
  return GEMGeometry_;
}

gemcrValidation::~gemcrValidation() {
}

auto_ptr<std::vector<TrajectorySeed> > gemcrValidation::findSeeds(MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits)
{
  auto_ptr<std::vector<TrajectorySeed> > tmptrajectorySeeds( new vector<TrajectorySeed>());
  for (auto hit1 : muRecHits){
    for (auto hit2 : muRecHits){
      if (hit1->globalPosition().y() < hit2->globalPosition().y()){
        LocalPoint segPos = hit1->localPosition();
        GlobalVector segDirGV(hit2->globalPosition().x() - hit1->globalPosition().x(),
                              (hit2->globalPosition().y() - hit1->globalPosition().y()),
                              hit2->globalPosition().z() - hit1->globalPosition().z());

        segDirGV *=10;
        LocalVector segDir = hit1->det()->toLocal(segDirGV);

        int charge= 1;
        LocalTrajectoryParameters param(segPos, segDir, charge);

        AlgebraicSymMatrix mat(5,0);
        mat = hit1->parametersError().similarityT( hit1->projectionMatrix() );
        LocalTrajectoryError error(asSMatrix<5>(mat));

        TrajectoryStateOnSurface tsos(param, error, hit1->det()->surface(), &*theService->magneticField());
        uint32_t id = hit1->rawId();
        PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState(tsos, id);

        edm::OwnVector<TrackingRecHit> seedHits;
        seedHits.push_back(hit1->hit()->clone());
        seedHits.push_back(hit2->hit()->clone());

        TrajectorySeed seed(seedTSOS,seedHits,alongMomentum);
        tmptrajectorySeeds->push_back(seed);
      }
    }
  }

  return tmptrajectorySeeds;
}

Trajectory gemcrValidation::makeTrajectory(TrajectorySeed seed, MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits, std::vector<GEMChamber> gemChambers, GEMChamber testChamber)
{
  PTrajectoryStateOnDet ptsd1(seed.startingState());
  DetId did(ptsd1.detId());
  const BoundPlane& bp = theService->trackingGeometry()->idToDet(did)->surface();
  TrajectoryStateOnSurface tsos = trajectoryStateTransform::transientState(ptsd1,&bp,&*theService->magneticField());
  TrajectoryStateOnSurface tsosCurrent = tsos;
  TransientTrackingRecHit::ConstRecHitContainer consRecHits;
  for (auto ch : gemChambers){
    if (ch == testChamber) continue;
    std::shared_ptr<MuonTransientTrackingRecHit> tmpRecHit;
    tsosCurrent = theService->propagator("SteppingHelixPropagatorAny")->propagate(tsosCurrent,ch.surface());
    if (!tsosCurrent.isValid()) continue;
    GlobalPoint tsosGP = tsosCurrent.freeTrajectoryState()->position();
    float maxR = 9999;
    for (auto hit : muRecHits){
      GEMDetId hitID(hit->rawId());
      if (hitID.chamberId() == ch.id() ){
        GlobalPoint hitGP = hit->globalPosition();
        if (abs(hitGP.x() - tsosGP.x()) > maxRes) continue;
        if (abs(hitGP.z() - tsosGP.z()) > 21.0) continue;
        float deltaR = (hitGP - tsosGP).mag();
        if (maxR > deltaR){
          tmpRecHit = hit;
          maxR = deltaR;
        }
      }
    }
    if (tmpRecHit){
      tsosCurrent = theUpdator->update(tsosCurrent, *tmpRecHit);
      consRecHits.push_back(tmpRecHit);
    }
  }
  if (consRecHits.size() <3) return Trajectory();
  vector<Trajectory> fitted = theSmoother->trajectories(seed, consRecHits, tsos);
  return fitted.front();
}

void gemcrValidation::analyze(const edm::Event& e, const edm::EventSetup& iSetup){

  theService->update(iSetup);

  edm::Handle<GEMRecHitCollection> gemRecHits;
  edm::Handle<vector<reco::Track>> gemTracks;
  edm::Handle<vector<TrajectorySeed>> trSeed;

  e.getByToken( this->InputTagToken_RH, gemRecHits);
  e.getByToken( this->InputTagToken_TR, gemTracks);
  e.getByToken( this->InputTagToken_TS, trSeed);
  if (!gemRecHits.isValid()) {
    edm::LogError("gemcrValidation") << "Cannot get strips by Token RecHits Token.\n";
    return ;
  }

  vector<bool> checkRH, checkTH, checkTR;
  vector<int> check_th_roll, check_th_vfat, check_tr_roll, check_tr_vfat;
  vector<double> tr_x, tr_y, th_x, th_y, x_err; 
  for(int c=0;c<n_ch;c++){
    checkRH.push_back(0);
    checkTH.push_back(0);
    checkTR.push_back(0);
    check_th_roll.push_back(-1);
    check_tr_roll.push_back(-1);
    check_th_vfat.push_back(-1);
    check_tr_vfat.push_back(-1);
    tr_x.push_back(-99);
    tr_y.push_back(-99);
    th_x.push_back(999);
    th_y.push_back(999);
    x_err.push_back(0.0);
  }
  for (GEMRecHitCollection::const_iterator recHit = gemRecHits->begin(); recHit != gemRecHits->end(); ++recHit){

    Float_t  rh_l_x = recHit->localPosition().x();
    Int_t  bx = recHit->BunchX();
    Int_t  clusterSize = recHit->clusterSize();
    Int_t  firstClusterStrip = recHit->firstClusterStrip();

    GEMDetId id((*recHit).gemId());
    int index = findIndex(id);
    checkRH[index] = 1;
    Short_t rh_roll = (Short_t) id.roll();
    LocalPoint recHitLP = recHit->localPosition();

    if ( GEMGeometry_->idToDet((*recHit).gemId()) == nullptr) {
      std::cout<<"This gem recHit did not matched with GEMGeometry."<<std::endl;
      continue;
    }
    GlobalPoint recHitGP = GEMGeometry_->idToDet((*recHit).gemId())->surface().toGlobal(recHitLP);
    Float_t     rh_g_X = recHitGP.x();
    Float_t     rh_g_Y = recHitGP.y();
    Float_t     rh_g_Z = recHitGP.z();
    int nVfat = 8*(findvfat(firstClusterStrip+clusterSize+1, 1, 128*3)-1) + (8-rh_roll);
    gem_chamber_x_y[index]->Fill(rh_l_x, rh_roll);
    gem_chamber_cl_size[index]->Fill(clusterSize, nVfat);
    gem_chamber_bx[index]->Fill(bx,rh_roll);
    gemcr_g->Fill(rh_g_X,rh_g_Z,rh_g_Y);
    gem_cls_tot->Fill(clusterSize);
    gem_bx_tot->Fill(bx);
    rh1_chamber->Fill(index);
    for(int i = firstClusterStrip; i < (firstClusterStrip + clusterSize); i++){
      gem_chamber_firedStrip[index]->Fill(i,rh_roll);
    }
    if (clusterSize < minCLS) continue;
    if (clusterSize > maxCLS) continue;
    rh2_chamber->Fill(index);
  }

  for (auto tch : gemChambers){
    MuonTransientTrackingRecHit::MuonRecHitContainer muRecHits;
    for (auto ch : gemChambers){
      if (tch == ch) continue;
      for (auto etaPart : ch.etaPartitions()){
        GEMDetId etaPartID = etaPart->id();
        GEMRecHitCollection::range range = gemRecHits->get(etaPartID);   
        for (GEMRecHitCollection::const_iterator rechit = range.first; rechit!=range.second; ++rechit){
          const GeomDet* geomDet(etaPart);
          if ((*rechit).clusterSize()<minCLS) continue;
          if ((*rechit).clusterSize()>maxCLS) continue;
          muRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));
        }
      }
    }
    auto_ptr<std::vector<TrajectorySeed> > trajectorySeeds( new vector<TrajectorySeed>());
    trajectorySeeds =findSeeds(muRecHits);
    Trajectory bestTrajectory;

    float maxChi2 = 1000;
    for (auto seed : *trajectorySeeds){
      Trajectory smoothed = makeTrajectory(seed, muRecHits, gemChambers,tch);
      if (smoothed.isValid()){
        //trajectorys->push_back(smoothed);
        if (maxChi2 > smoothed.chiSquared()/float(smoothed.ndof())){
          maxChi2 = smoothed.chiSquared()/float(smoothed.ndof());
          bestTrajectory = smoothed;
        }
      }
    }





  }
  if (gemTracks.isValid()){
    tr_size->Fill(gemTracks.product()->size());
    for(auto tr: *gemTracks.product()){
      for(trackingRecHit_iterator hit = tr.recHitsBegin(); hit != tr.recHitsEnd(); hit++){
        LocalPoint tlp = (*hit)->localPosition();
        GEMDetId tid = GEMDetId((*hit)->geographicalId());
        int index = findIndex(tid);
        //GlobalPoint recHitGP = gemChambers[index].surface().toGlobal(tlp);
        checkTR[index] = 1;
        tr_x[index] = tlp.x();
        tr_y[index] = tid.roll();
        check_tr_roll[index] = tid.roll();
        float min_x = gemChambers[index].etaPartition(tid.roll())->centreOfStrip(0).x();
        float max_x = gemChambers[index].etaPartition(tid.roll())->centreOfStrip(128*3).x();
        check_tr_vfat[index] = findvfat(tlp.x(),min_x, max_x);
        x_err[index] = (*hit)->localPositionError().xx();
        rh3_chamber->Fill(index);
        //gem_chamber_residual[index]->Fill(tlp.x(), tid.roll());
        //cout << (*hit)->clusterSize() << endl;
        //cout << "roll : " << tid.roll() << ", vfat : " << findvfat(tr_rh_lp.x(),min_x, max_x) << endl;
      }
      TrajectorySeed seed = trSeed.product()->back();
      PTrajectoryStateOnDet ptsd1(seed.startingState());
      DetId did(ptsd1.detId());
      const BoundPlane& bp = theService->trackingGeometry()->idToDet(did)->surface();
      TrajectoryStateOnSurface tsos = trajectoryStateTransform::transientState(ptsd1,&bp,&*theService->magneticField());
      TrajectoryStateOnSurface tsosCurrent = tsos;
      for(int c=0; c<n_ch;c++){
        GEMChamber ch = gemChambers[c];
        tsosCurrent = theService->propagator("SteppingHelixPropagatorAny")->propagate(tsosCurrent,ch.surface());
        if (!tsosCurrent.isValid()) continue;
        Global3DPoint gtrp = tsosCurrent.freeTrajectoryState()->position();
        Local3DPoint tlp = ch.surface().toLocal(gtrp);
        if (!ch.surface().bounds().inside(tlp)){continue;}
        int n_roll = ch.nEtaPartitions();
        double minDely = 50.;
        int mRoll = -1;
        for (int r=0;r<n_roll;r++){
          Local3DPoint rtlp = ch.etaPartition(r+1)->surface().toLocal(gtrp);
          if(minDely > abs(rtlp.y())){minDely = abs(rtlp.y()); mRoll = r+1;}
        }
        if (mRoll == -1){cout << "no mRoll" << endl;continue;}
        int n_strip = ch.etaPartition(mRoll)->nstrips();
        double min_x = ch.etaPartition(mRoll)->centreOfStrip(0).x();
        double max_x = ch.etaPartition(mRoll)->centreOfStrip(n_strip).x();
        if ((tlp.x()>(min_x)) & (tlp.x() < (max_x)) ){ 
          check_th_roll[c] = mRoll;
          check_th_vfat[c] = findvfat(tlp.x(), min_x, max_x);
          checkTH[c] = 1; 
          th_x[c] = tlp.x();
          if (checkTR[c]) { 
            Local3DPoint rtlp = ch.etaPartition(tr_y[c])->surface().toLocal(gtrp);
            th_y[c] = rtlp.y();
          }
          continue;
        }
      }
      /*int flag = 0;
      for(int c=0;c<n_ch;c++){
        if (checkTR[c] and checkTH[c]){
          if (abs(th_x[c]-tr_x[c])<1.0) {flag++;}
          else {checkTR[c] = 0; checkTH[c] = 0;  }
        }
      }
      if (flag<3){ continue; }*/
      /*vector<int> workingCh = {2,3,4,5};
      bool flag1 = 1;
      for (auto x : workingCh){
        if (!checkTH[x]) flag1 = 0;
      }
      if (!flag1) continue;
      cout << "all TrackHit" << endl;
      bool flagTest = 1;
      for (auto x : workingCh){ 
        cout << x << " : " << checkTR[x] << endl;
        if (!checkTR[x]) flagTest = 0;
      }
      if (flagTest) cout << "all trRecHit" << endl;
      for (auto x : workingCh){
        bool flag2 = 1;
        for (auto y : workingCh){
          if (x != y) {
            if (!checkTR[y]) flag2 = 0;
          }
        }
        if (!flag2) continue;
        gem_chamber_th2D_eff[x]->Fill(th_x[x], check_th_roll[x]);
        if (checkTR[x]) {gem_chamber_tr2D_eff[x]->Fill(th_x[x], check_th_roll[x]);}
                     
      }*/

      for (int tc=0; tc<n_ch;tc++){
        if (checkTR[tc]) {
              gem_chamber_local_x[tc]->Fill(tr_x[tc],th_x[tc]);
              gem_chamber_residual[tc]->Fill(th_x[tc], tr_x[tc]-th_x[tc]);
        }
      }
      for(int c=0;c<n_ch;c++){
        //if (!checkTH[c] and checkTR[c]) {cout << "bad track" << endl;}
        if (!checkTH[c]) {continue;}
        if (!checkTR[c]){
          gem_chamber_th2D_eff[c]->Fill(check_th_vfat[c],check_th_roll[c]);
          gem_chamber_thxy_eff[c]->Fill(th_x[c],check_th_roll[c]);
        }
        if(checkTR[c]){
          gem_chamber_th2D_eff[c]->Fill(check_th_vfat[c],check_th_roll[c]);
          gem_chamber_thxy_eff[c]->Fill(th_x[c],check_th_roll[c]);

          gem_chamber_tr2D_eff[c]->Fill(check_th_vfat[c],check_th_roll[c]);
          gem_chamber_trxy_eff[c]->Fill(th_x[c],check_th_roll[c]);
          //gem_chamber_residual[c]->Fill(tr_x[c]-th_x[c], th_y[c]);
          //gem_chamber_residual_r[c]->Fill( sqrt((tr_x[c]-th_x[c])*(tr_x[c]-th_x[c])+ (tr_y[c]-th_y[c])*(tr_y[c]-th_y[c])));
          //gem_chamber_residual_r[c]->Fill( th_x[c] - tr_x[c]);
          //gem_chamber_residual_r[c]->Fill( sqrt((th_x[c]-tr_x[c])*(th_x[c]-tr_x[c]) + th_y[c]*th_y[c]));
          //gem_chamber_residual_r[c]->Fill(th_x[c]-tr_x[c]);
        }
      }
    }
  }
}


