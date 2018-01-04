/** \class GEMCosmicMuonForQC8
 * Produces a collection of tracks's in GEM cosmic ray stand. 
 *
 * \author Jason Lee
 */

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

using namespace std;

class GEMCosmicMuonForQC8 : public edm::stream::EDProducer<> {
public:
  /// Constructor
  explicit GEMCosmicMuonForQC8(const edm::ParameterSet&);
  /// Destructor
  virtual ~GEMCosmicMuonForQC8() {}
  /// Produce the GEMSegment collection
  void produce(edm::Event&, const edm::EventSetup&) override;
  double maxCLS;
  double minCLS;
  double trackChi2, trackResX, trackResY;
  double MulSigmaOnWindow;
private:
  int iev; // events through
  edm::EDGetTokenT<GEMRecHitCollection> theGEMRecHitToken;
  CosmicMuonSmoother* theSmoother;
  MuonServiceProxy* theService;
  KFUpdator* theUpdator;
  unique_ptr<std::vector<TrajectorySeed> > findSeeds(MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits);
  int findSeeds(std::vector<TrajectorySeed> *tmptrajectorySeeds, 
    MuonTransientTrackingRecHit::MuonRecHitContainer &seedupRecHits, 
    MuonTransientTrackingRecHit::MuonRecHitContainer &seeddnRecHits, 
    std::vector<int> &vecnTypeHitsUp, std::vector<int> &vecnTypeHitsDn);
  Trajectory makeTrajectory(TrajectorySeed seed, MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits, std::vector<GEMChamber> gemChambers, GEMChamber testChamber);
};

GEMCosmicMuonForQC8::GEMCosmicMuonForQC8(const edm::ParameterSet& ps) : iev(0) {
  maxCLS = ps.getParameter<double>("maxClusterSize");
  minCLS = ps.getParameter<double>("minClusterSize");
  trackChi2 = ps.getParameter<double>("trackChi2");
  trackResX = ps.getParameter<double>("trackResX");
  trackResY = ps.getParameter<double>("trackResY");
  MulSigmaOnWindow = ps.getParameter<double>("MulSigmaOnWindow");
  theGEMRecHitToken = consumes<GEMRecHitCollection>(ps.getParameter<edm::InputTag>("gemRecHitLabel"));
  // register what this produces
  edm::ParameterSet serviceParameters = ps.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  edm::ParameterSet smootherPSet = ps.getParameter<edm::ParameterSet>("MuonSmootherParameters");
  theSmoother = new CosmicMuonSmoother(smootherPSet,theService);
  theUpdator = new KFUpdator();
  produces<reco::TrackCollection>();
  produces<TrackingRecHitCollection>();
  produces<reco::TrackExtraCollection>();
  produces<vector<Trajectory> >();
  produces<vector<TrajectorySeed> >();
  produces<vector<int> >();
  produces<vector<unsigned int> >();
}

void GEMCosmicMuonForQC8::produce(edm::Event& ev, const edm::EventSetup& setup) {
  //  cout << "GEMCosmicMuonForQC8::start producing segments for " << ++iev << "th event with gem data" << endl;
  
  unique_ptr<reco::TrackCollection >          trackCollection( new reco::TrackCollection() );
  unique_ptr<TrackingRecHitCollection >       trackingRecHitCollection( new TrackingRecHitCollection() );
  unique_ptr<reco::TrackExtraCollection >     trackExtraCollection( new reco::TrackExtraCollection() );
  unique_ptr<vector<Trajectory> >             trajectorys( new vector<Trajectory>() );
  unique_ptr<vector<TrajectorySeed> >         trajectorySeeds( new vector<TrajectorySeed>() );
  unique_ptr<vector<int> >                    trajectoryChIdx( new vector<int>() );
  unique_ptr<vector<unsigned int> >           trajectoryType( new vector<unsigned int>() );
  unique_ptr<vector<double> >                 trajectorySeedsHits( new vector<double>() );
  TrackingRecHitRef::key_type recHitsIndex = 0;
  TrackingRecHitRefProd recHitCollectionRefProd = ev.getRefBeforePut<TrackingRecHitCollection>();
  reco::TrackExtraRef::key_type trackExtraIndex = 0;
  reco::TrackExtraRefProd trackExtraCollectionRefProd = ev.getRefBeforePut<reco::TrackExtraCollection>();
  
  theService->update(setup);

  edm::ESHandle<GEMGeometry> gemg;
  setup.get<MuonGeometryRecord>().get(gemg);
  const GEMGeometry* mgeom = &*gemg;
  
  vector<GEMChamber> gemChambers;
  //int n_ch;
  
  const std::vector<const GEMSuperChamber*>& superChambers_ = mgeom->superChambers();   
  for (auto sch : superChambers_){
    int n_lay = sch->nChambers();
    for (int l=0;l<n_lay;l++){
      gemChambers.push_back(*sch->chamber(l+1));
    }
  }
  //n_ch = gemChambers.size();

  //vector<const DetLayer*> gemLayers = theMuonLayers->allDTLayers();
  
  // get the collection of GEMRecHit
  edm::Handle<GEMRecHitCollection> gemRecHits;
  ev.getByToken(theGEMRecHitToken,gemRecHits);

  if (gemRecHits->size() <3){
    ev.put(std::move(trajectorySeeds));
    ev.put(std::move(trackCollection));
    ev.put(std::move(trackingRecHitCollection));
    ev.put(std::move(trackExtraCollection));
    ev.put(std::move(trajectorys));
    ev.put(std::move(trajectoryChIdx));
    ev.put(std::move(trajectoryType));
    return;
  }
  
  //  cout << "GEMCosmicMuonForQC8::gemRecHits " << gemRecHits->size() << endl;
  
  vector<double> vecChamType = {1,3, 0,0, 0,0, 0,0, 4,2, 
                                1,3, 0,0, 0,0, 0,0, 4,2, 
                                1,3, 0,0, 0,0, 0,0, 4,2};
  
  int countTC = 0;
    
  for (auto tch : gemChambers){
    countTC++;
    MuonTransientTrackingRecHit::MuonRecHitContainer muRecHits;
    MuonTransientTrackingRecHit::MuonRecHitContainer seedupRecHits;
    MuonTransientTrackingRecHit::MuonRecHitContainer seeddnRecHits;
    
    std::vector<int> vecnTypeHitsUp;
    std::vector<int> vecnTypeHitsDn;
    
    int nUpType = 2;
    int nDnType = 1;
    
    int nIdxTestCh = tch.id().chamber() + tch.id().layer() - 2;
    //int nColTestCh = nIdxTestCh / 10;
    int nIsTopBottom = 0;
    
    //if ( vecChamType[ nIdxTestCh ] == 2 ) {nUpType = 4; nIsTopBottom = 1;}
    //if ( vecChamType[ nIdxTestCh ] == 1 ) {nDnType = 3; nIsTopBottom = 1;}
    if ( vecChamType[ nIdxTestCh ] == 2 || vecChamType[ nIdxTestCh ] == 1 ) {
      nUpType = 4;
      nDnType = 3;
      nIsTopBottom = 1;
    }
    
    int countC = 0;
    int TCN = 0;
    for (auto ch : gemChambers){
      countC += 1;
      if (tch == ch) continue;
      int nHitOnceFilter = 0;
      for (auto etaPart : ch.etaPartitions()){
        GEMDetId etaPartID = etaPart->id();
        GEMRecHitCollection::range range = gemRecHits->get(etaPartID);   
        int nRoll = (int)etaPartID.roll();
        for (GEMRecHitCollection::const_iterator rechit = range.first; rechit!=range.second; ++rechit){
          const GeomDet* geomDet(etaPart);
          if ((*rechit).clusterSize()<minCLS) continue;
          if ((*rechit).clusterSize()>maxCLS) continue;
          muRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));
          
          if ( nHitOnceFilter == 0 ) {
            TCN++;
            nHitOnceFilter = 1;
          }
          
          int nIdxCh = ch.id().chamber() + ch.id().layer() - 2;
          //int nColCh = nIdxCh / 10;
          if ( vecChamType[ nIdxCh ] == nUpType ) {
            seedupRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));
            if ( nIsTopBottom != 0 && nRoll == 1/* && nColCh == nColTestCh*/ ) {
              vecnTypeHitsUp.push_back(1);
            } else if ( nIsTopBottom != 0 && nRoll == 8/* && nColCh == nColTestCh*/ ) {
              vecnTypeHitsUp.push_back(8);
            } else {
              vecnTypeHitsUp.push_back(0);
            }
          } else if ( vecChamType[ nIdxCh ] == nDnType ) {
            seeddnRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));
            if ( nIsTopBottom != 0 && nRoll == 1/* && nColCh == nColTestCh*/ ) {
              vecnTypeHitsDn.push_back(1);
            } else if ( nIsTopBottom != 0 && nRoll == 8/* && nColCh == nColTestCh*/ ) {
              vecnTypeHitsDn.push_back(8);
            } else {
              vecnTypeHitsDn.push_back(0);
            }
          }
        }
      }
    }
    if (muRecHits.size()<3){ 
      //cout << "tracking passed due to # recHit under 3" << endl;
      //if ( testRecHits.size() > 0 ) gem_chamber_track[findIndex(tch.id())]->Fill(6.5);
      continue;}
    if (TCN < 3) continue;
    //if ( testRecHits.size() > 0 ) gem_chamber_track[findIndex(tch.id())]->Fill(0.5);
    vector<TrajectorySeed> trajSeedsBody;
    std::vector<TrajectorySeed> *trajSeeds = &trajSeedsBody;
    int nNumNormalSeed;
    nNumNormalSeed = findSeeds(trajSeeds, seedupRecHits, seeddnRecHits, vecnTypeHitsUp, vecnTypeHitsDn);
    Trajectory bestTrajectory;
    TrajectorySeed bestSeed;
    //trajectoryh->Fill(0, trajSeeds->size());
    float maxChi2 = 10000000.0;
    int countTR = 0;
    int nIdxBest = -1;
    for (auto seed : *trajSeeds){
      Trajectory smoothed = makeTrajectory(seed, muRecHits, gemChambers,tch);
      
      countTR += 1;
      if (smoothed.isValid()){
        //trajectoryh->Fill(2,1);
        float dProbChiNDF = smoothed.chiSquared()/float(smoothed.ndof());
        //float dProbChiNDF = TMath::Prob(smoothed.chiSquared(), smoothed.ndof());
        if (maxChi2 > dProbChiNDF){
          maxChi2 = dProbChiNDF;
          bestTrajectory = smoothed;
          bestSeed = seed;
          nIdxBest = countTR - 1;
        }
      }else{
        //strKeep += TString::Format("non valid\n");
        //trajectoryh->Fill(1,1);
        //cout << "trajectory " << countTR << " is not valid" << endl;
      }
    }
    
    if (!bestTrajectory.isValid()) continue; //{cout<<"no Best Trajectory" << endl; continue;}
    //if ( testRecHits.size() > 0 ) gem_chamber_track[findIndex(tch.id())]->Fill(1.5);
    //if (maxChi2 > trackChi2) continue;
    if (maxChi2 > trackChi2) {
      continue;
    }
    
    unsigned int unTypeSeed = 0;
    
    if ( nIdxBest >= nNumNormalSeed ) {
      unTypeSeed = 1;
    }
    
    const FreeTrajectoryState* ftsAtVtx = bestTrajectory.geometricalInnermostState().freeState();
    
    GlobalPoint pca = ftsAtVtx->position();
    math::XYZPoint persistentPCA(pca.x(),pca.y(),pca.z());
    GlobalVector p = ftsAtVtx->momentum();
    math::XYZVector persistentMomentum(p.x(),p.y(),p.z());
    
    reco::Track track(bestTrajectory.chiSquared(), 
                      bestTrajectory.ndof(true),
                      persistentPCA,
                      persistentMomentum,
                      ftsAtVtx->charge(),
                      ftsAtVtx->curvilinearError());
   
    // create empty collection of Segments
    //cout << "GEMCosmicMuon::track " << track.pt() << endl;
    
    // reco::TrackExtra tx(track.outerPosition(), track.outerMomentum(), track.outerOk(),
    // 		      track.innerPosition(), track.innerMomentum(), track.innerOk(),
    // 		      track.outerStateCovariance(), track.outerDetId(),
    // 		      track.innerStateCovariance(), track.innerDetId(),
    // 		      track.seedDirection(), edm::RefToBase<TrajectorySeed>());
    // 		      //, bestTrajectory.seedRef() );
    reco::TrackExtra tx;
    //tx.setResiduals(track.residuals());
    //adding rec hits
    Trajectory::RecHitContainer transHits = bestTrajectory.recHits();
    unsigned int nHitsAdded = 0;
    for (Trajectory::RecHitContainer::const_iterator recHit = transHits.begin(); recHit != transHits.end(); ++recHit) {
      TrackingRecHit *singleHit = (**recHit).hit()->clone();
      trackingRecHitCollection->push_back( singleHit );  
      ++nHitsAdded;
    }
    tx.setHits(recHitCollectionRefProd, recHitsIndex, nHitsAdded);
    recHitsIndex +=nHitsAdded;
    
    trackExtraCollection->push_back(tx );
    
    reco::TrackExtraRef trackExtraRef(trackExtraCollectionRefProd, trackExtraIndex++ );
    track.setExtra(trackExtraRef);
    trackCollection->push_back(track);
    
    //trajectoryh->Fill(3,1);
    //gem_chamber_bestChi2[findIndex(tch.id())]->Fill(maxChi2);
    //if ( testRecHits.size() > 0 ) gem_chamber_track[findIndex(tch.id())]->Fill(2.5);
    trajectorys->push_back(bestTrajectory);
    trajectorySeeds->push_back(bestSeed);
    trajectoryChIdx->push_back(countTC - 1);
    trajectoryType->push_back(unTypeSeed);
  }
  
  // fill the collection
  // put collection in event
  ev.put(std::move(trajectorySeeds));
  ev.put(std::move(trackCollection));
  ev.put(std::move(trackingRecHitCollection));
  ev.put(std::move(trackExtraCollection));
  ev.put(std::move(trajectorys));
  ev.put(std::move(trajectoryChIdx));
  ev.put(std::move(trajectoryType));
  
}

unique_ptr<std::vector<TrajectorySeed> > GEMCosmicMuonForQC8::findSeeds(MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits)
{
  unique_ptr<std::vector<TrajectorySeed> > trajectorySeeds( new vector<TrajectorySeed>());
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
        //float p_err = 0.2;
        //mat[0][0]= p_err;
        LocalTrajectoryError error(asSMatrix<5>(mat));

        // get first hit
        TrajectoryStateOnSurface tsos(param, error, hit1->det()->surface(), &*theService->magneticField());
        //cout << "GEMCosmicMuonForQC8::tsos " << tsos << endl;
        uint32_t id = hit1->rawId();
           PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState(tsos, id);
        
        edm::OwnVector<TrackingRecHit> seedHits;
        seedHits.push_back(hit1->hit()->clone());
        seedHits.push_back(hit2->hit()->clone());

        TrajectorySeed seed(seedTSOS,seedHits,alongMomentum);
        trajectorySeeds->push_back(seed);
      }
    }
  }
  
  return trajectorySeeds;
}


int GEMCosmicMuonForQC8::findSeeds(std::vector<TrajectorySeed> *tmptrajectorySeeds, 
    MuonTransientTrackingRecHit::MuonRecHitContainer &seedupRecHits, 
    MuonTransientTrackingRecHit::MuonRecHitContainer &seeddnRecHits, 
    std::vector<int> &vecnTypeHitsUp, std::vector<int> &vecnTypeHitsDn)
{
  std::vector<TrajectorySeed> tmptrajectorySeedsEdgeEta;
  
  int nIdxHitUp, nIdxHitDn;
  
  nIdxHitDn = -1;
  
  for (auto hit1 : seeddnRecHits){
    nIdxHitDn++;
    nIdxHitUp = -1;
    for (auto hit2 : seedupRecHits){
      nIdxHitUp++;
      if (hit1->globalPosition().y() < hit2->globalPosition().y()) {
        LocalPoint segPos = hit1->localPosition();
        GlobalVector segDirGV(hit2->globalPosition().x() - hit1->globalPosition().x(),
                              (hit2->globalPosition().y() - hit1->globalPosition().y()),
                              hit2->globalPosition().z() - hit1->globalPosition().z());

        segDirGV *=10;
        //segDirGV *=1;
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
        
        if ( ( vecnTypeHitsUp[ nIdxHitUp ] == 1 && vecnTypeHitsDn[ nIdxHitDn ] == 1 ) || 
            ( vecnTypeHitsUp[ nIdxHitUp ] == 8 && vecnTypeHitsDn[ nIdxHitDn ] == 8 ) )
        {
          tmptrajectorySeedsEdgeEta.push_back(seed);
        } else {
          tmptrajectorySeeds->push_back(seed);
        }
      }
    }
  }
  
  int nSizeNormalSeed = (int)tmptrajectorySeeds->size();
  
  for ( auto seedEdgeEta : tmptrajectorySeedsEdgeEta ) {
    tmptrajectorySeeds->push_back(seedEdgeEta);
  }
  return nSizeNormalSeed;
}


/*float gemcrValidation::CalcWindowWidthX(GPSeed *pVecSeed, GlobalPoint *pPCurr) {
  if ( 1 != 1 ) {
    return 5.0;
  }
  
  float fDev = trackResX * MulSigmaOnWindow;
  
  float fYCenterSeed = 0.5 * ( pVecSeed->P1.y() + pVecSeed->P2.y() );
  float fYDiffSeed = pVecSeed->P2.y() - fYCenterSeed;
  float fYDiffHit = pPCurr->y() - fYCenterSeed;
  
  //printf("%lf\n", abs(fYDiffHit / fYDiffSeed));
  return max(fDev * abs(fYDiffHit / fYDiffSeed), fDev);
}


float gemcrValidation::CalcWindowWidthY(GPSeed *pVecSeed, GlobalPoint *pPCurr) { // Y : local
  if ( 1 != 1 ) {
    return 5.0;
  }
  
  float fDev = trackResY * MulSigmaOnWindow;
  
  float fYCenterSeed = 0.5 * ( pVecSeed->P1.y() + pVecSeed->P2.y() );
  float fYDiffSeed = pVecSeed->P2.y() - fYCenterSeed;
  float fYDiffHit = pPCurr->y() - fYCenterSeed;
  
  //printf("%lf\n", abs(fYDiffHit / fYDiffSeed));
  return max(fDev * abs(fYDiffHit / fYDiffSeed), fDev);
}*/


Trajectory GEMCosmicMuonForQC8::makeTrajectory(TrajectorySeed seed, MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits, std::vector<GEMChamber> gemChambers, GEMChamber testChamber)
{
  PTrajectoryStateOnDet ptsd1(seed.startingState());
  DetId did(ptsd1.detId());
  const BoundPlane& bp = theService->trackingGeometry()->idToDet(did)->surface();
  TrajectoryStateOnSurface tsos = trajectoryStateTransform::transientState(ptsd1,&bp,&*theService->magneticField());
  TrajectoryStateOnSurface tsosCurrent = tsos;
  TransientTrackingRecHit::ConstRecHitContainer consRecHits;
  for (auto ch : gemChambers){
    //if (ch == testChamber) continue;
    std::shared_ptr<MuonTransientTrackingRecHit> tmpRecHit;
    tsosCurrent = theService->propagator("SteppingHelixPropagatorAny")->propagate(tsosCurrent,ch.surface());
    if (!tsosCurrent.isValid()) return Trajectory();
    GlobalPoint tsosGP = tsosCurrent.freeTrajectoryState()->position();
    float maxR = 9999;
    for (auto hit : muRecHits){
      GEMDetId hitID(hit->rawId());
      if (hitID.chamberId() == ch.id() ){
        GlobalPoint hitGP = hit->globalPosition();

        //double y_err = hit->localPositionError().yy();
        //double x_err = hit->localPositionError().xx();
        if (abs(hitGP.x() - tsosGP.x()) > trackResX * MulSigmaOnWindow) continue;
        if (abs(hitGP.z() - tsosGP.z()) > trackResY * MulSigmaOnWindow * hit->localPositionError().yy() ) continue; // global z, local y
        //if (abs(hitGP.z() - tsosGP.z()) > y_err*trackResY) continue;
        float deltaR = (hitGP - tsosGP).mag();
        if (maxR > deltaR){
          tmpRecHit = hit;
          maxR = deltaR;
        }
      }
    }
    if (tmpRecHit){
      //tsosCurrent = theUpdator->update(tsosCurrent, *tmpRecHit);
      consRecHits.push_back(tmpRecHit);
      //GlobalPoint recHitGP = GEMGeometry_->idToDet((*tmpRecHit).rawId())->surface().toGlobal(tmpRecHit->localPosition());
    }
  }
  if (consRecHits.size() <3) return Trajectory();
  vector<Trajectory> fitted = theSmoother->trajectories(seed, consRecHits, tsos);
  if ( fitted.size() <= 0 ) return Trajectory();
  return fitted.front();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GEMCosmicMuonForQC8);
