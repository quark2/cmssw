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
  MuonServiceProxy* theService;
  
};

GEMCosmicMuon::GEMCosmicMuon(const edm::ParameterSet& ps) : iev(0) {
  theGEMRecHitToken = consumes<GEMRecHitCollection>(ps.getParameter<edm::InputTag>("gemRecHitLabel"));
  // register what this produces
  edm::ParameterSet serviceParameters = ps.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  edm::ParameterSet smootherPSet = ps.getParameter<edm::ParameterSet>("MuonSmootherParameters");
  theSmoother = new CosmicMuonSmoother(smootherPSet,theService);
  produces<vector<reco::Track> >();
  produces<vector<Trajectory> >();
}

void GEMCosmicMuon::produce(edm::Event& ev, const edm::EventSetup& setup) {
  cout << "GEMCosmicMuon::start producing segments for " << ++iev << "th event with gem data" << endl;
  std::auto_ptr<vector<reco::Track> > oc( new vector<reco::Track> );
  std::auto_ptr<vector<Trajectory> > ocTra( new vector<Trajectory> );

  theService->update(setup);

  edm::ESHandle<MagneticField> theField;  
  setup.get<IdealMagneticFieldRecord>().get(theField);

  edm::ESHandle<GEMGeometry> gemg;
  setup.get<MuonGeometryRecord>().get(gemg);
  const GEMGeometry* mgeom = &*gemg;

  //vector<const DetLayer*> gemLayers = theMuonLayers->allDTLayers();
  
  // get the collection of GEMRecHit
  edm::Handle<GEMRecHitCollection> gemRecHits;
  ev.getByToken(theGEMRecHitToken,gemRecHits);

  if (gemRecHits->size() <3){
     ev.put(oc);
     return;
  }
    
  cout << "GEMCosmicMuon::gemRecHits " << gemRecHits->size() << endl;
  
  MuonTransientTrackingRecHit::MuonRecHitContainer muRecHits;
  edm::OwnVector<TrackingRecHit> ovTrackRecHits;
  TransientTrackingRecHit::ConstRecHitContainer consRecHits;

  //for (auto ch : chambers){
  // Create the chamber Id
  for (int i =0; i < 10; i++){
    for (int j=1; j < 3; j++){
      GEMDetId chamberID(1,1,1, j, i, 0);
      const GEMChamber* ch = mgeom->chamber(chamberID);
      if (!ch) continue;
      for (auto etaPart : ch->etaPartitions()){
	GEMDetId etaPartID = etaPart->id();
	// Get the GEM-Segment which relies on this chamber
	GEMRecHitCollection::range range = gemRecHits->get(etaPartID);

	//cout<< "Number of GEM rechits available , from chamber: "<< etaPartID<<endl;
	// Create the MuonTransientTrackingRecHit
	for (GEMRecHitCollection::const_iterator rechit = range.first; rechit!=range.second; ++rechit){
	  const GeomDet* geomDet(etaPart);	
	  muRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));

	  //GEMRecHit *newRH = rechit->clone();
	  //TransientTrackingRecHit::ConstRecHitPointer trkRecHit(newRH);
	  consRecHits.push_back(muRecHits.back());
	  ovTrackRecHits.push_back(rechit->clone());
	  cout << "GEMCosmicMuon::hit " << muRecHits.back()->globalPosition() << endl;
	}
      }
    }
  }
  cout<< "Number of GEM rechits   = " << muRecHits.size()<<endl;
  cout<< "Number of GEM container = " << ovTrackRecHits.size()<<endl;
  
  // cout << "GEMCosmicMuon::gemRecHits " << gemRecHits->size() << endl;
  // GEMRecHit* gemHit = 0;
  // for (GEMRecHitCollection::const_iterator recHit = gemRecHits->begin(); recHit != gemRecHits->end(); ++recHit){
  //   container.push_back(recHit->clone());
  //   gemHit = recHit->clone();
  // }
  // //if (!gemHit) return;
  
  if (ovTrackRecHits.size() < 3){
    ev.put(ocTra);
    ev.put(oc);
    return;
  }
  
  MuonTransientTrackingRecHit::MuonRecHitPointer hit = muRecHits[0];
  MuonTransientTrackingRecHit::MuonRecHitPointer hit2 = muRecHits[1];

  // cout << "GEMCosmicMuon::hit " << gemHit->gemId() << endl;
  // cout << "GEMCosmicMuon::firstClusterStrip " << gemHit->firstClusterStrip() << endl;
  // cout << "GEMCosmicMuon::clusterSize " << gemHit->clusterSize() << endl;
  cout << "GEMCosmicMuon::hit " << hit->localPosition() << endl;
  cout << "GEMCosmicMuon::hit " << hit->globalPosition() << endl;
  
  LocalPoint segPos = hit->localPosition();
  GlobalVector segDirGV(hit2->globalPosition().x() - hit->globalPosition().x(),
			(hit2->globalPosition().y() - hit->globalPosition().y()),
			hit2->globalPosition().z() - hit->globalPosition().z());

  segDirGV *=10;
  LocalVector segDir = hit->det()->toLocal(segDirGV);
  cout << "GEMCosmicMuon::GlobalVector " << segDirGV << endl;
  cout << "GEMCosmicMuon::LocalVector  " << segDir << endl;
  
  int charge= 1;
  LocalTrajectoryParameters param(segPos, segDir, charge);
  
  AlgebraicSymMatrix mat(5,0);
  mat = hit->parametersError().similarityT( hit->projectionMatrix() );
  float p_err = 0.2;
  mat[0][0]= p_err;
  LocalTrajectoryError error(asSMatrix<5>(mat));

  // get first hit
  TrajectoryStateOnSurface tsos(param, error, hit->det()->surface(), &*theField);
  //cout << "GEMCosmicMuon::tsos " << tsos << endl;
  uint32_t id = hit->rawId();
   
  PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState(tsos, id);
  TrajectorySeed seed(seedTSOS,ovTrackRecHits,alongMomentum);
  
  vector<Trajectory> fitted = theSmoother->trajectories(seed, consRecHits, tsos);

  if (fitted.empty()) return;
  Trajectory smoothed = fitted.front();
  cout << "GEMCosmicMuon::Trajectory " << smoothed.foundHits() << endl;
  
  if (smoothed.foundHits() < 3) return;
  
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
  cout << "GEMCosmicMuon::track " << track.pt() << endl;


  ocTra->push_back(smoothed);
  oc->push_back(track);
  // fill the collection
  // put collection in event
  ev.put(ocTra);
  ev.put(oc);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GEMCosmicMuon);
