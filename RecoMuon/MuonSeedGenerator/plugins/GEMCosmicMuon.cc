/** \class GEMCosmicMuon derived by CSCSegmentProducer 
 * Produces a collection of GEMSegment's in endcap muon GEMs. 
 *
 * \author Piet Verwilligen
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

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/CosmicMuonProducer/interface/CosmicMuonSmoother.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

using namespace std;

class GEMCosmicMuon : public edm::stream::EDProducer<> {
public:
  /// Constructor
  explicit GEMCosmicMuon(const edm::ParameterSet&);
  /// Destructor
  virtual ~GEMCosmicMuon() {}
  /// Produce the GEMSegment collection
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  int iev; // events through
  edm::EDGetTokenT<GEMRecHitCollection> theGEMRecHitToken;
  CosmicMuonSmoother* theSmoother;
  const MuonServiceProxy* theService;
  
};

GEMCosmicMuon::GEMCosmicMuon(const edm::ParameterSet& ps) : iev(0) {
  theGEMRecHitToken = consumes<GEMRecHitCollection>(ps.getParameter<edm::InputTag>("gemRecHitLabel"));
  // register what this produces
  edm::ParameterSet serviceParameters = ps.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  edm::ParameterSet smootherPSet = ps.getParameter<edm::ParameterSet>("MuonSmootherParameters");
  theSmoother = new CosmicMuonSmoother(smootherPSet,theService);
  produces<vector<reco::Track> >();
}

void GEMCosmicMuon::produce(edm::Event& ev, const edm::EventSetup& setup) {
  LogDebug("GEMCosmicMuon") << "start producing segments for " << ++iev << "th event with gem data";
  // find the geometry (& conditions?) for this event & cache it in the builder
  // edm::ESHandle<GEMGeometry> gemg;
  // setup.get<MuonGeometryRecord>().get(gemg);
  // const GEMGeometry* mgeom = &*gemg;
  
  edm::ESHandle<MagneticField> theField;  
  setup.get<IdealMagneticFieldRecord>().get(theField);

  // get the collection of GEMRecHit
  edm::Handle<GEMRecHitCollection> gemRecHits;
  ev.getByToken(theGEMRecHitToken,gemRecHits);
  
  edm::OwnVector<TrackingRecHit> container;
  GEMRecHit* gemHit = 0;
  for (GEMRecHitCollection::const_iterator recHit = gemRecHits->begin(); recHit != gemRecHits->end(); ++recHit){
    container.push_back(recHit->clone());
    gemHit = recHit->clone();
  }
  if (!gemHit) return;
  
  MuonTransientTrackingRecHit::MuonRecHitPointer hit = MuonTransientTrackingRecHit::specificBuild( gemHit->det(), gemHit->hit() );

  LocalPoint segPos = hit->localPosition();
  LocalVector segDir(0,0,-1.);
  int charge= 0;
  LocalTrajectoryParameters param(segPos,segDir, charge);

  AlgebraicSymMatrix mat(5,0);
  mat = hit->parametersError().similarityT( hit->projectionMatrix() );
  float p_err = 0.2;
  mat[0][0]= p_err;
  LocalTrajectoryError error(asSMatrix<5>(mat));
  
  // get first hit
  TrajectoryStateOnSurface tsos(param, error, hit->det()->surface(), &*theField);
  uint32_t id = hit->rawId();
    
  PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState(tsos, id);
  TrajectorySeed seed(seedTSOS,container,alongMomentum);

  const Trajectory* theTraj = new Trajectory(seed,alongMomentum);

  Trajectory smoothed = theSmoother->trajectory(*theTraj); 


  const FreeTrajectoryState* ftsAtVtx = smoothed.geometricalInnermostState().freeState();

  GlobalPoint pca = ftsAtVtx->position();
  math::XYZPoint persistentPCA(pca.x(),pca.y(),pca.z());
  GlobalVector p = ftsAtVtx->momentum();
  math::XYZVector persistentMomentum(p.x(),p.y(),p.z());
  
  reco::Track track(smoothed.chiSquared(), 
		    smoothed.ndof(true),
		    persistentPCA,
		    persistentMomentum,
		    ftsAtVtx->charge(),
		    ftsAtVtx->curvilinearError());
 
  // create empty collection of Segments
  std::auto_ptr<vector<reco::Track> > oc( new vector<reco::Track> );
  oc->push_back(track);
  // fill the collection
  // put collection in event
  ev.put(oc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GEMCosmicMuon);
