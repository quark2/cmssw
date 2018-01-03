#include "Validation/GEMCR/interface/gemcrValidation.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/GEMGeometry/interface/GEMSuperChamber.h"
#include <Geometry/GEMGeometry/interface/GEMGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>

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

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"

//#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include <iomanip>

#include <TCanvas.h>
using namespace std;
gemcrValidation::gemcrValidation(const edm::ParameterSet& cfg): GEMBaseValidation(cfg)
{
  isMC = cfg.getParameter<bool>("isMC");
  InputTagToken_RH = consumes<GEMRecHitCollection>(cfg.getParameter<edm::InputTag>("recHitsInputLabel"));
  InputTagToken_TR = consumes<vector<reco::Track>>(cfg.getParameter<edm::InputTag>("tracksInputLabel"));
  InputTagToken_TS = consumes<vector<TrajectorySeed>>(cfg.getParameter<edm::InputTag>("seedInputLabel"));
  InputTagToken_TI = consumes<vector<int>>(cfg.getParameter<edm::InputTag>("chNoInputLabel"));
  InputTagToken_TT = consumes<vector<unsigned int>>(cfg.getParameter<edm::InputTag>("seedTypeInputLabel"));
  InputTagToken_DG = consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("gemDigiLabel"));
  //InputTagToken_ST = consumes<vector<SimTrack>>(cfg.getParameter<edm::InputTag>("simTrack"));
  if ( isMC ) InputTagToken_US = consumes<edm::HepMCProduct>(cfg.getParameter<edm::InputTag>("genVtx"));
  edm::ParameterSet serviceParameters = cfg.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  minCLS = cfg.getParameter<double>("minClusterSize"); 
  maxCLS = cfg.getParameter<double>("maxClusterSize");
  maxRes = cfg.getParameter<double>("maxResidual");
  makeTrack = cfg.getParameter<bool>("makeTrack");
  trackChi2 = cfg.getParameter<double>("trackChi2");
  trackResY = cfg.getParameter<double>("trackResY"); 
  trackResX = cfg.getParameter<double>("trackResX");
  MulSigmaOnWindow = cfg.getParameter<double>("MulSigmaOnWindow");
  
  if ( isMC ) {
    /*fScinHPosZ   = cfg.getParameter<double>("ScintilUpperZ");
    fScinHLeft   = cfg.getParameter<double>("ScintilUpperLeft");
    fScinHRight  = cfg.getParameter<double>("ScintilUpperRight");
    fScinHTop    = cfg.getParameter<double>("ScintilUpperTop");
    fScinHBottom = cfg.getParameter<double>("ScintilUpperBottom");
    
    fScinLPosZ   = cfg.getParameter<double>("ScintilLowerZ");
    fScinLLeft   = cfg.getParameter<double>("ScintilLowerLeft");
    fScinLRight  = cfg.getParameter<double>("ScintilLowerRight");
    fScinLTop    = cfg.getParameter<double>("ScintilLowerTop");
    fScinLBottom = cfg.getParameter<double>("ScintilLowerBottom");*/
    
    ScinUPosZ.push_back(cfg.getParameter<double>("ScintilUpper00Z"));
    ScinUXMin.push_back(cfg.getParameter<double>("ScintilUpper00XMin"));
    ScinUXMax.push_back(cfg.getParameter<double>("ScintilUpper00XMax"));
    ScinUYMin.push_back(cfg.getParameter<double>("ScintilUpper00YMin"));
    ScinUYMax.push_back(cfg.getParameter<double>("ScintilUpper00YMax"));
    
    ScinLPosZ.push_back(cfg.getParameter<double>("ScintilLower00Z"));
    ScinLXMin.push_back(cfg.getParameter<double>("ScintilLower00XMin"));
    ScinLXMax.push_back(cfg.getParameter<double>("ScintilLower00XMax"));
    ScinLYMin.push_back(cfg.getParameter<double>("ScintilLower00YMin"));
    ScinLYMax.push_back(cfg.getParameter<double>("ScintilLower00YMax"));
    
    ScinUPosZ.push_back(cfg.getParameter<double>("ScintilUpper01Z"));
    ScinUXMin.push_back(cfg.getParameter<double>("ScintilUpper01XMin"));
    ScinUXMax.push_back(cfg.getParameter<double>("ScintilUpper01XMax"));
    ScinUYMin.push_back(cfg.getParameter<double>("ScintilUpper01YMin"));
    ScinUYMax.push_back(cfg.getParameter<double>("ScintilUpper01YMax"));
    
    ScinLPosZ.push_back(cfg.getParameter<double>("ScintilLower01Z"));
    ScinLXMin.push_back(cfg.getParameter<double>("ScintilLower01XMin"));
    ScinLXMax.push_back(cfg.getParameter<double>("ScintilLower01XMax"));
    ScinLYMin.push_back(cfg.getParameter<double>("ScintilLower01YMin"));
    ScinLYMax.push_back(cfg.getParameter<double>("ScintilLower01YMax"));
    
    ScinUPosZ.push_back(cfg.getParameter<double>("ScintilUpper02Z"));
    ScinUXMin.push_back(cfg.getParameter<double>("ScintilUpper02XMin"));
    ScinUXMax.push_back(cfg.getParameter<double>("ScintilUpper02XMax"));
    ScinUYMin.push_back(cfg.getParameter<double>("ScintilUpper02YMin"));
    ScinUYMax.push_back(cfg.getParameter<double>("ScintilUpper02YMax"));
    
    ScinLPosZ.push_back(cfg.getParameter<double>("ScintilLower02Z"));
    ScinLXMin.push_back(cfg.getParameter<double>("ScintilLower02XMin"));
    ScinLXMax.push_back(cfg.getParameter<double>("ScintilLower02XMax"));
    ScinLYMin.push_back(cfg.getParameter<double>("ScintilLower02YMin"));
    ScinLYMax.push_back(cfg.getParameter<double>("ScintilLower02YMax"));
    
    ScinUPosZ.push_back(cfg.getParameter<double>("ScintilUpper03Z"));
    ScinUXMin.push_back(cfg.getParameter<double>("ScintilUpper03XMin"));
    ScinUXMax.push_back(cfg.getParameter<double>("ScintilUpper03XMax"));
    ScinUYMin.push_back(cfg.getParameter<double>("ScintilUpper03YMin"));
    ScinUYMax.push_back(cfg.getParameter<double>("ScintilUpper03YMax"));
    
    ScinLPosZ.push_back(cfg.getParameter<double>("ScintilLower03Z"));
    ScinLXMin.push_back(cfg.getParameter<double>("ScintilLower03XMin"));
    ScinLXMax.push_back(cfg.getParameter<double>("ScintilLower03XMax"));
    ScinLYMin.push_back(cfg.getParameter<double>("ScintilLower03YMin"));
    ScinLYMax.push_back(cfg.getParameter<double>("ScintilLower03YMax"));
    
    ScinUPosZ.push_back(cfg.getParameter<double>("ScintilUpper04Z"));
    ScinUXMin.push_back(cfg.getParameter<double>("ScintilUpper04XMin"));
    ScinUXMax.push_back(cfg.getParameter<double>("ScintilUpper04XMax"));
    ScinUYMin.push_back(cfg.getParameter<double>("ScintilUpper04YMin"));
    ScinUYMax.push_back(cfg.getParameter<double>("ScintilUpper04YMax"));
    
    ScinLPosZ.push_back(cfg.getParameter<double>("ScintilLower04Z"));
    ScinLXMin.push_back(cfg.getParameter<double>("ScintilLower04XMin"));
    ScinLXMax.push_back(cfg.getParameter<double>("ScintilLower04XMax"));
    ScinLYMin.push_back(cfg.getParameter<double>("ScintilLower04YMin"));
    ScinLYMax.push_back(cfg.getParameter<double>("ScintilLower04YMax"));
  }
  
  edm::ParameterSet smootherPSet = cfg.getParameter<edm::ParameterSet>("MuonSmootherParameters");
  theSmoother = new CosmicMuonSmoother(smootherPSet, theService);
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

  //gemcr_g = ibooker.book3D("gemcr_g","GEMCR GLOBAL RECHITS", 200,-100,100,20,-55,55, 140,0,140);
  gemcr_g = ibooker.book3D("gemcr_g","GEMCR GLOBAL RECHITS", 260,-130,130,30,-82.5,82.5, 140,0,140);
  gemcrTr_g = ibooker.book3D("gemcrTr_g","GEMCR GLOBAL RECHITS", 260,-130,130,30,-82.5,82.5, 140,0,140);
  gemcrCf_g = ibooker.book3D("gemcrCf_g","GEMCR GLOBAL RECHITS CONFIRMED", 260,-130,130,30,-82.5,82.5, 140,0,140);
  gemcrTrScint_g = ibooker.book3D("gemcrTrScint_g","GEMCR GLOBAL RECHITS", 260,-130,130,30,-82.5,82.5, 140,0,140);
  gemcrCfScint_g = ibooker.book3D("gemcrCfScint_g","GEMCR GLOBAL RECHITS SCINTILLATED", 260,-130,130,30,-82.5,82.5, 140,0,140);
  gem_cls_tot = ibooker.book1D("cluster_size","Cluseter size",20,0,20);  
  gem_bx_tot = ibooker.book1D("bx", "BX" , 30, -15,15);
  tr_size = ibooker.book1D("tr_size", "track size",10,0,10);
  tr_hit_size = ibooker.book1D("tr_hit_size", "hit size in track",15,0,15); 
  trajectoryh = ibooker.book1D("trajectory","trajectory", 4,0,4);
  trajectoryh->setBinLabel(1, "total seeds");
  trajectoryh->setBinLabel(2, "unvalid");
  trajectoryh->setBinLabel(3, "valid");
  trajectoryh->setBinLabel(4, "passed chi2");
  firedMul = ibooker.book1D("firedMul","fired chamber multiplicity",n_ch+1,0,n_ch+1);
  firedChamber = ibooker.book1D("firedChamber", "fired chamber",n_ch,0,n_ch);
  scinUpperHit = ibooker.book3D("scinUpperHit","UPPER SCINTILLATOR GLOBAL", 260,-130,130,30,-82.5,82.5, 140,0,140);
  scinLowerHit = ibooker.book3D("scinUpperHit","LOWER SCINTILLATOR GLOBAL", 260,-130,130,30,-82.5,82.5, 140,0,140);
  scinUpperRecHit = ibooker.book3D("scinLowerRecHit","UPPER SCINTILLATOR GLOBAL RECHITS", 260,-130,130,30,-82.5,82.5, 140,0,140);
  scinLowerRecHit = ibooker.book3D("scinLowerRecHit","LOWER SCINTILLATOR GLOBAL RECHITS", 260,-130,130,30,-82.5,82.5, 140,0,140);
  
  resXSim = ibooker.book1D("residualx_sim", " residual x (sim)",200,-3,3);
  resYByErrSim = ibooker.book1D("residualy_sim", " residual y (sim)",200,-3,3);
  hitXErr = ibooker.book1D("x_err", "x_err",200,0,5);
  hitYErr = ibooker.book1D("y_err", "y_err",200,0,50);

  tr_chamber = ibooker.book1D("tr_eff_ch", "tr rec /chamber",n_ch,1,n_ch); 
  th_chamber = ibooker.book1D("th_eff_ch", "tr hit/chamber",n_ch,0,n_ch); 
  rh_chamber = ibooker.book1D("rh_eff_ch", "rec hit/chamber",n_ch,0,n_ch); 
  rh1_chamber = ibooker.book1D("rh1_chamber", "all recHits",n_ch,0,n_ch); 
  rh2_chamber = ibooker.book1D("rh2_chamber", "cut passed recHits ",n_ch,0,n_ch); 
  rh3_chamber = ibooker.book1D("rh3_chamber", "tracking recHits",n_ch,0,n_ch); 
  
  rh3_chamber_scint = ibooker.book1D("rh3_chamber_scint", "tracking recHits, scintillated",n_ch,0,n_ch); 
  
  events_withtraj = ibooker.book1D("events_withtraj", "# of events with a reconstructed trajectory",8,0,8); 
  
  events_withtraj->setBinLabel(1, "events with firedCh >= 4");
  events_withtraj->setBinLabel(2, "events with firedCh >= 5");
  events_withtraj->setBinLabel(3, "events with firedCh >= 6");
  events_withtraj->setBinLabel(4, "events with top-bottom seed");
  events_withtraj->setBinLabel(5, "events with tr, firedCh >= 4");
  events_withtraj->setBinLabel(6, "events with tr, firedCh >= 5");
  events_withtraj->setBinLabel(7, "events with tr, firedCh >= 6");
  events_withtraj->setBinLabel(8, "events with tr, top-bottom seed");
  
  aftershots = ibooker.book2D("plain_aftershots", "aftershots", 25, -130, 130, 25, -82.5, 82.5); 
  
  //diffTrajGenRec = ibooker.book3D("diffTrajGenRec","Difference between GEN track and RECO track", 
  //  100, 0, 3, 100, -1, 1, 100, 0, 50);
  
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
   firedChamber->setBinLabel(c+1,b_name);
   rh3_chamber_scint->setBinLabel(c+1,b_name);
  }
  for(int c = 0; c<n_ch;c++){
     GEMDetId gid = gemChambers[c].id();
     string h_name = "chamber_"+to_string(gid.chamber())+"_layer_"+to_string(gid.layer());
     gem_chamber_x_y.push_back(ibooker.book2D(h_name+"_recHit",h_name+" recHit", 500,-25,25,8,1,9));
     gem_chamber_cl_size.push_back(ibooker.book2D(h_name+"_recHit_size", h_name+" recHit size", 50,0,50,24,0,24));
     gem_chamber_bx.push_back(ibooker.book2D(h_name+"_bx", h_name+" BX", 30,-15,15,10,0,10));
     gem_chamber_pad_vfat.push_back(ibooker.book2D(h_name+"_pad_vfat", h_name+" pad (vfat)", 3,1,4,8,1,9));
     gem_chamber_copad_vfat.push_back(ibooker.book2D(h_name+"_copad_vfat", h_name+" copad (vfat)", 3,1,4,8,1,9));
     gem_chamber_pad_vfat_withmul.push_back(ibooker.book2D(h_name+"_pad_vfat_withmul", h_name+" pad (vfat), with multiplicity", 3,1,4,8,1,9));
     gem_chamber_copad_vfat_withmul.push_back(ibooker.book2D(h_name+"_copad_vfat_withmul", h_name+" copad (vfat), with multiplicity", 3,1,4,8,1,9));
     gem_chamber_tr2D_eff.push_back(ibooker.book2D(h_name+"_recHit_efficiency", h_name+" recHit efficiency", 3,1,4,8,1,9));
     gem_chamber_th2D_eff.push_back(ibooker.book2D(h_name+"_th2D_eff", h_name+"_th2D_eff", 3,1,4,8,1,9));
     gem_chamber_trxroll_eff.push_back(ibooker.book2D(h_name+"_trxroll_eff", h_name+" recHit efficiency", 50,-25,25,8,1,9));
     gem_chamber_trxy_eff.push_back(ibooker.book2D(h_name+"_trxy_eff", h_name+" recHit efficiency", 50,-25,25,120,-60,60));
     gem_chamber_thxroll_eff.push_back(ibooker.book2D(h_name+"_thxroll_eff", h_name+"_th2D_eff", 50,-25,25,8,1,9));
     gem_chamber_thxy_eff.push_back(ibooker.book2D(h_name+"_thxy_eff", h_name+"_th2D_eff", 50,-25,25,120,-60,60));
     gem_chamber_residual.push_back(ibooker.book2D(h_name+"_residual", h_name+" residual", 500,-25,25,100,-5,5));
     gem_chamber_residualX1DSim.push_back(ibooker.book1D(h_name+"_residualX1DSim", h_name+" residual x 1D (sim)", 200,-3,3));
     gem_chamber_residualY1DSim.push_back(ibooker.book1D(h_name+"_residualY1DSim", h_name+" residual y 1D (sim)", 200,-3,3));
     gem_chamber_local_x.push_back(ibooker.book2D(h_name+"_local_x", h_name+" local x",500,-25,25,500,-25,25));
     gem_chamber_digi_digi.push_back(ibooker.book2D(h_name+"_digi_gemDigi", h_name+" gemDigi (DIGI)", 384,0,384,8,1,9));
     gem_chamber_digi_recHit.push_back(ibooker.book2D(h_name+"_recHit_gemDigi", h_name+" gemDigi (recHit)", 384,0,384,8,1,9));
     gem_chamber_digi_CLS.push_back(ibooker.book2D(h_name+"_CLS_gemDigi", h_name+" gemDigi (CLS)", 384,0,384,8,1,9));
     gem_chamber_hitMul.push_back(ibooker.book1D(h_name+"_hit_mul", h_name+" hit multiplicity",25,0,25 ));
     gem_chamber_vfatHitMul.push_back(ibooker.book2D(h_name+"_vfatHit_mul", h_name+" vfat hit multiplicity",25,0,25, 24,0,24));
     gem_chamber_stripHitMul.push_back(ibooker.book2D(h_name+"_stripHit_mul", h_name+" strip hit multiplicity", 150,0,150,9,0,9));
     gem_chamber_bestChi2.push_back(ibooker.book1D(h_name+"_bestChi2", h_name+" #chi^{2} distribution", trackChi2*10,0,trackChi2));
     gem_chamber_track.push_back(ibooker.book1D(h_name+"_track", h_name+" track",7,0,7));
     
     gem_chamber_th2D_eff_scint.push_back(ibooker.book2D(h_name+"_th2D_eff_scint", h_name+"_th2D_eff, scintillated", 3,1,4,8,1,9));
     gem_chamber_thxroll_eff_scint.push_back(ibooker.book2D(h_name+"_thxroll_eff_scint", h_name+"_th2D_eff, scintillated", 50,-25,25,8,1,9));
     gem_chamber_thxy_eff_scint.push_back(ibooker.book2D(h_name+"_thxy_eff_scint", h_name+"_th2D_eff, scintillated", 50,-25,25,120,-60,60));
     
     gem_chamber_tr2D_eff_scint.push_back(ibooker.book2D(h_name+"_recHit_efficiency_scint", h_name+" recHit efficiency, scintillated", 3,1,4,8,1,9));
     gem_chamber_trxroll_eff_scint.push_back(ibooker.book2D(h_name+"_trxroll_eff_scint", h_name+" recHit efficiency, scintillated", 50,-25,25,8,1,9));
     gem_chamber_trxy_eff_scint.push_back(ibooker.book2D(h_name+"_trxy_eff_scint", h_name+" recHit efficiency, scintillated", 50,-25,25,120,-60,60));
     gem_chamber_local_x_scint.push_back(ibooker.book2D(h_name+"_local_x_scint", h_name+" local x, scintillated",500,-25,25,500,-25,25));
     gem_chamber_residual_scint.push_back(ibooker.book2D(h_name+"_residual_scint", h_name+" residual, scintillated", 500,-25,25,100,-5,5));
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

////////////////////////////////////////////////////////////////////////////////
// To turn on or turn off methods for bad trajectory veto, seek them: 
// CONE_WINDOW
// VETO_TOO_FEW_RECHITS
// VETO_TOO_SHORT_SEED
////////////////////////////////////////////////////////////////////////////////

int g_nEvt = 0;
int g_nNumRecHit = 0;
int g_nNumFiredCh = 0;
int g_nNumFiredChValid = 0;
int g_nNumTrajHit = 0;
int g_nNumTrajHit6 = 0;
int g_nNumMatched = 0;

double g_dMinY = 100000.0, g_dMaxY = -10000000.0;

gemcrValidation::~gemcrValidation() {
  printf("res1 : %i\n", g_nEvt);
  printf("res2 : %i\n", g_nNumRecHit);
  printf("res3 : %i\n", g_nNumFiredCh);
  printf("res4 : %i\n", g_nNumFiredChValid);
  printf("res5 : %i\n", g_nNumTrajHit6);
  printf("res6 : %i\n", g_nNumTrajHit);
  printf("res7 : %i\n", g_nNumMatched);
  
  printf("Y range : %lf, %lf\n", g_dMinY, g_dMaxY);
}


bool gemcrValidation::isPassedScintillators(GlobalPoint trajGP1, GlobalPoint trajGP2) {
  unsigned int i;
  
  // Getting the direction of trajectory by using the hits in the seed
  double dTrajDirX = trajGP2.x() - trajGP1.x();
  double dTrajDirY = trajGP2.y() - trajGP1.y();
  double dTrajDirZ = trajGP2.z() - trajGP1.z();
  
  int nIsHitUpper, nIsHitLower;
  
  double dXUHit = 0.0, dYUHit = 0.0, dZUHit = 0.0;
  double dXLHit = 0.0, dYLHit = 0.0, dZLHit = 0.0;
  
  if ( -0.0000001 < dTrajDirY && dTrajDirY < 0.0000001 ) return false; // for safety
  
  nIsHitUpper = 0;
  nIsHitLower = 0;
  
  for ( i = 0 ; i < ScinUPosZ.size() ; i++ ) {
    // Time to solve 1-order equation
    double dTTraj = ( ScinUPosZ[ i ] - trajGP1.y() ) / dTrajDirY;
    
    // Finding the coordinates of the hit position on the y-plane containing the upper scintillator
    dXUHit = trajGP1.x() + dTTraj * dTrajDirX;
    dYUHit = ScinUPosZ[ i ];
    dZUHit = trajGP1.z() + dTTraj * dTrajDirZ;
    
    // The above lines are probably not needed to be calculated in each iteration
    // because all entries of ScinUPosZ may be same.
    // But for... safety...
    
    // Checking whether the hit position in the upper scintillator
    if ( ( ScinUXMin[ i ] <= dXUHit && dXUHit <= ScinUXMax[ i ]  ) && 
       ( ScinUYMin[ i ]  <= dZUHit && dZUHit <= ScinUYMax[ i ] ) )
    {
      nIsHitUpper = 1;
      break;
    }
  }
  
  for ( i = 0 ; i < ScinLPosZ.size() ; i++ ) {
    // Time to solve 1-order equation
    double dTTraj = ( ScinLPosZ[ i ] - trajGP1.y() ) / dTrajDirY;
    
    // Finding the coordinates of the hit position on the y-plane containing the upper scintillator
    dXLHit = trajGP1.x() + dTTraj * dTrajDirX;
    dYLHit = ScinLPosZ[ i ];
    dZLHit = trajGP1.z() + dTTraj * dTrajDirZ;
    
    // The above lines are probably not needed to be calculated in each iteration
    // because all entries of ScinLPosZ may be same.
    // But for... safety...
    
    // Checking whether the hit position in the upper scintillator
    if ( ( ScinLXMin[ i ] <= dXLHit && dXLHit <= ScinLXMax[ i ]  ) && 
       ( ScinLYMin[ i ]  <= dZLHit && dZLHit <= ScinLYMax[ i ] ) )
    {
      nIsHitLower = 1;
      break;
    }
  }
  
  /*vecPosHit.push_back(dXUHit);
  vecPosHit.push_back(dYUHit);
  vecPosHit.push_back(dZUHit);
  
  vecPosHit.push_back(dXLHit);
  vecPosHit.push_back(dYLHit);
  vecPosHit.push_back(dZLHit);*/
  
  scinUpperHit->Fill(dXUHit, dZUHit, dYUHit);
  scinLowerHit->Fill(dXLHit, dZLHit, dYLHit);
  
  if ( nIsHitUpper == 1 ) scinUpperRecHit->Fill(dXUHit, dZUHit, dYUHit);
  if ( nIsHitLower == 1 ) scinLowerRecHit->Fill(dXLHit, dZLHit, dYLHit);
  
  return ( nIsHitUpper == 1 && nIsHitLower == 1 );
}


int g_nNumTest = 0;


void gemcrValidation::analyze(const edm::Event& e, const edm::EventSetup& iSetup){
  g_nEvt++;

  //std::cout << "analyze on!" << std::endl;
  theService->update(iSetup);

  edm::Handle<GEMRecHitCollection> gemRecHits;
  if (!isMC){
    edm::Handle<GEMDigiCollection> digis;
    e.getByToken( this->InputTagToken_DG, digis);
    int nNumDigi = 0;
    for (GEMDigiCollection::DigiRangeIterator gemdgIt = digis->begin(); gemdgIt != digis->end(); ++gemdgIt){
      nNumDigi++;
      const GEMDetId& gemId = (*gemdgIt).first;
      int index = findIndex(gemId);
      const GEMDigiCollection::Range& range = (*gemdgIt).second;
      for ( auto digi = range.first; digi != range.second; ++digi ) {
        gem_chamber_digi_digi[index]->Fill(digi->strip(),gemId.roll());
      }
    }
    //std::cout << "num of digi : " << nNumDigi << std::endl;
  }
  e.getByToken( this->InputTagToken_RH, gemRecHits);
  if (!gemRecHits.isValid()) {
    edm::LogError("gemcrValidation") << "Cannot get strips by Token RecHits Token.\n";
    return ;
  }
  
  float fXGenGP1x = 0.0, fXGenGP1y = 0.0, fXGenGP1z = 0.0;
  float fXGenGP2x = 0.0, fXGenGP2y = 0.0, fXGenGP2z = 0.0;
  
  int nNumCurrFiredCh = 0;
  
  if ( isMC ) {
    edm::Handle<edm::HepMCProduct> genVtx;
    e.getByToken( this->InputTagToken_US, genVtx);
    HepMC::GenParticle *genMuon = genVtx->GetEvent()->barcode_to_particle(1);
    
    double dUnitGen = 0.1;
    
    fXGenGP1x = dUnitGen * genMuon->production_vertex()->position().x();
    fXGenGP1y = dUnitGen * genMuon->production_vertex()->position().y();
    fXGenGP1z = dUnitGen * genMuon->production_vertex()->position().z();
    
    fXGenGP2x = fXGenGP1x + dUnitGen * genMuon->momentum().x();
    fXGenGP2y = fXGenGP1y + dUnitGen * genMuon->momentum().y();
    fXGenGP2z = fXGenGP1z + dUnitGen * genMuon->momentum().z();
    
    //printf("%lf, %lf, %lf ---> %lf, %lf, %lf\n", fXGenGP1x, fXGenGP1z, fXGenGP1y, fXGenGP2x, fXGenGP2z, fXGenGP2y);
    
    Float_t fVecX/*, fVecY*/, fVecZ;
    int arrnFired[ 32 ] = {0, };
    
    fVecX = genMuon->momentum().x() / genMuon->momentum().y();
    //fVecY = 1.0;
    fVecZ = genMuon->momentum().z() / genMuon->momentum().y();
    
    for ( GEMRecHitCollection::const_iterator recHit = gemRecHits->begin(); recHit != gemRecHits->end(); ++recHit ) {
      int nIdxCh = findIndex((*recHit).gemId());
      GlobalPoint recHitGP = GEMGeometry_->idToDet((*recHit).gemId())->surface().toGlobal(recHit->localPosition());
      
      if ( arrnFired[ nIdxCh ] == 0 ) {
        arrnFired[ nIdxCh ] = 1;
        //g_nNumFiredCh++;
        nNumCurrFiredCh++;
      }
      
      Float_t fDiffY = recHitGP.y() - fXGenGP1y;
      
      Float_t fXGenHitX = fXGenGP1x + fDiffY * fVecX;
      Float_t fXGenHitZ = fXGenGP1z + fDiffY * fVecZ;
      //Float_t fXGenHitZ = recHitGP.y();
      
      //printf("  recHit : %i, (%0.5f, %0.5f, %0.5f) <... (%0.5f, %0.5f, %0.5f)\n", nIdxCh, 
      //    recHitGP.x(), recHitGP.z(), recHitGP.y(), fXGenHitX, fXGenHitZ, recHitGP.y());
      
      //resXSim->Fill(recHitGP.x() - fXGenHitX);
      resXSim->Fill(recHitGP.x() - fXGenHitX);
      gem_chamber_residualX1DSim[ nIdxCh ]->Fill(recHitGP.x() - fXGenHitX);
      hitXErr->Fill(recHit->localPositionError().xx());
      
      resYByErrSim->Fill(( recHitGP.z() - fXGenHitZ ) / recHit->localPositionError().yy());
      gem_chamber_residualY1DSim[ nIdxCh ]->Fill(( recHitGP.z() - fXGenHitZ ) / recHit->localPositionError().yy());
      hitYErr->Fill(recHit->localPositionError().yy());
      
      g_nNumRecHit++;
    }
    
    aftershots->Fill(fXGenGP1x + 400.0 * fVecX, fXGenGP1z + 400.0 * fVecZ);
    
    g_nNumFiredCh += nNumCurrFiredCh;
    if ( nNumCurrFiredCh > 6 ) g_nNumFiredChValid += nNumCurrFiredCh;
  }
  
  GlobalPoint genGPos1(fXGenGP1x, fXGenGP1y, fXGenGP1z);
  GlobalPoint genGPos2(fXGenGP2x, fXGenGP2y, fXGenGP2z);
  
  vector<bool> firedCh;
  vector<int> rMul;
  vector<vector<int>> vMul(n_ch, vector<int>(24, 0));
  vector<vector<int>> sMul(n_ch, vector<int>(9, 0));
  for (int c=0;c<n_ch;c++){
    firedCh.push_back(0);
    rMul.push_back(0);
  }
  //std::cout << "num of hit : " << gemRecHits->size() << std::endl;
  TString strListRecHit("");
  for (GEMRecHitCollection::const_iterator recHit = gemRecHits->begin(); recHit != gemRecHits->end(); ++recHit){

    Float_t  rh_l_x = recHit->localPosition().x();
    Int_t  bx = recHit->BunchX();
    Int_t  clusterSize = recHit->clusterSize();
    Int_t  firstClusterStrip = recHit->firstClusterStrip();

    GEMDetId id((*recHit).gemId());
    int index = findIndex(id);
    firedCh[index] = 1;
    rMul[index] += 1;
    //checkRH[index] = 1;
    Short_t rh_roll = (Short_t) id.roll();
    LocalPoint recHitLP = recHit->localPosition();
    //printf("%i, %lf, %i, %i\n", index, rh_l_x, clusterSize, firstClusterStrip);

    if ( GEMGeometry_->idToDet((*recHit).gemId()) == nullptr) {
      std::cout<<"This gem recHit did not matched with GEMGeometry."<<std::endl;
      continue;
    }
    GlobalPoint recHitGP = GEMGeometry_->idToDet((*recHit).gemId())->surface().toGlobal(recHitLP);
    Float_t     rh_g_X = recHitGP.x();
    Float_t     rh_g_Y = recHitGP.y();
    Float_t     rh_g_Z = recHitGP.z();
    int nVfat = 8*(findvfat(firstClusterStrip+clusterSize*0.5, 0, 128*3)-1) + (8-rh_roll);
    vMul[index][nVfat] += 1;
    gem_chamber_x_y[index]->Fill(rh_l_x, rh_roll);
    gem_chamber_cl_size[index]->Fill(clusterSize, nVfat);
    gem_chamber_bx[index]->Fill(bx,rh_roll);
    gemcr_g->Fill(rh_g_X,rh_g_Z,rh_g_Y);
    gem_cls_tot->Fill(clusterSize);
    gem_bx_tot->Fill(bx);
    rh1_chamber->Fill(index);
    for(int i = firstClusterStrip; i < (firstClusterStrip + clusterSize); i++){
      gem_chamber_digi_recHit[index]->Fill(i,rh_roll);
    }
    strListRecHit += TString::Format("recHit : %lf, %lf, %lf (%i, %i) (%s)\n", rh_g_X, rh_g_Z, rh_g_Y, index + 1 - ( index % 2 ), ( index % 2 ) + 1, ( clusterSize >= minCLS && clusterSize <= maxCLS ? "okay" : "cut" ));
    if (clusterSize < minCLS) continue;
    if (clusterSize > maxCLS) continue;
    rh2_chamber->Fill(index);
    if ( rMul[ index ] <= 1 ) gem_chamber_track[ index ]->Fill(3.5);
    for(int i = firstClusterStrip; i < (firstClusterStrip + clusterSize); i++){
      sMul[index][rh_roll] +=1;
      sMul[index][0] +=1;
      gem_chamber_digi_CLS[index]->Fill(i,rh_roll);
    }
    //printf("recHit : %lf, %lf, %lf (%i, %i)\n", rh_g_X, rh_g_Z, rh_g_Y, index + 1 - ( index % 2 ), ( index % 2 ) + 1);
  }
  
  for ( int ich = 0 ; ich < n_ch ; ich++ ) {
    for ( int ivfat = 0 ; ivfat < 24 ; ivfat++ ) {
      if ( vMul[ ich ][ ivfat ] > 0 ) {
        int nRoll = 8 - ivfat % 8;
        int nVFat = ivfat / 8 + 1;
        
        gem_chamber_pad_vfat[ ich ]->Fill(nVFat, nRoll);
        for ( int mulidx = 0 ; mulidx < vMul[ ich ][ ivfat ] ; mulidx++ ) gem_chamber_pad_vfat_withmul[ ich ]->Fill(nVFat, nRoll);
        
        if ( vMul[ ich ^ 0x1 ][ ivfat ] > 0 ) {
          gem_chamber_copad_vfat[ ich ]->Fill(nVFat, nRoll);
          for ( int mulidx = 0 ; mulidx < vMul[ ich ^ 0x1 ][ ivfat ] ; mulidx++ ) gem_chamber_copad_vfat_withmul[ ich ]->Fill(nVFat, nRoll);
        }
      }
    }
  }
  
  int fChMul = 0;
  for(int c=0;c<n_ch;c++){
    gem_chamber_hitMul[c]->Fill(rMul[c]);
    for(int v=0; v<24;v++){
      gem_chamber_vfatHitMul[c]->Fill(vMul[c][v],v);    
    }
    for(int r=0; r<9;r++){
      gem_chamber_stripHitMul[c]->Fill(sMul[c][r],r);
    } 
    if (firedCh[c]){ 
      firedChamber->Fill(c+0.5);
      fChMul += 1;
    }
  }
  if (fChMul <= 3) //cout << "less then 3 chambers fired !"<< endl;
  {
    cout << "less then 3 chambers fired !"<< endl;
    /*for (GEMRecHitCollection::const_iterator recHit = gemRecHits->begin(); recHit != gemRecHits->end(); ++recHit){
      GlobalPoint recHitGP = GEMGeometry_->idToDet((*recHit).gemId())->surface().toGlobal(recHit->localPosition());
      printf("  recHit : (%0.5f, %0.5f, %0.5f)\n", 
        recHitGP.x(), recHitGP.z(), recHitGP.y());
    }*/
  }
  firedMul->Fill(fChMul);
    
  bool bIsScint = isPassedScintillators(genGPos1, genGPos2);
  
  /// Tracking start
  vector<double> chamSetEff = {0,0, 0,0, 0,0, 0,0, 0,0, 
                              //97.0,95.0, 96.0,97.0, 97.0,96.0, 98.0,96.0, 96.0,97.0, 
                              97.0,97.0, 97.0,97.0, 97.0,97.0, 97.0,97.0, 97.0,97.0, 
                              //97.0,97.0, 97.0,97.0, 97.0,97.0, 97.0,97.0, 0,0, 
                              //97.0,97.0, 97.0,97.0, 97.0,97.0, 97.0,97.0, 97.0,97.0, 
                              0,0, 0,0, 0,0, 0,0, 0,0};
  
  vector<double> vecChamType = {1,3, 0,0, 0,0, 0,0, 4,2, 
                                1,3, 0,0, 0,0, 0,0, 4,2, 
                                1,3, 0,0, 0,0, 0,0, 4,2};
  
  int fCha = 10;
  int lCha = 19;
  
  edm::Handle<std::vector<int>> idxChTraj;
  e.getByToken( this->InputTagToken_TI, idxChTraj);
  
  edm::Handle<std::vector<TrajectorySeed>> seedGCM;
  e.getByToken( this->InputTagToken_TS, seedGCM);
  
  edm::Handle<std::vector<unsigned int>> seedTypes;
  e.getByToken( this->InputTagToken_TT, seedTypes);
  
  if ( idxChTraj->size() == 0 ) return;

  if (!makeTrack) return; 
  int countTC = 0;
  int nIsTraceGet = 0;
  int nIsLongSeed = 0;
  //TString strKeep;
  if ( fChMul == 4 ) events_withtraj->Fill(0.5);
  if ( fChMul == 5 ) events_withtraj->Fill(1.5);
  if ( fChMul >= 6 ) events_withtraj->Fill(2.5);
  for (auto tch : gemChambers){
    //strKeep = TString("");
    countTC += 1;
    MuonTransientTrackingRecHit::MuonRecHitContainer testRecHits;
    if (isMC){if (tch.id().chamber()<9 and tch.id().chamber()>20) continue;}
    for (auto etaPart : tch.etaPartitions()){
      GEMDetId etaPartID = etaPart->id();
      GEMRecHitCollection::range range = gemRecHits->get(etaPartID);
      for (GEMRecHitCollection::const_iterator rechit = range.first; rechit!=range.second; ++rechit){
          const GeomDet* geomDet(etaPart);
          if ((*rechit).clusterSize()<minCLS) continue;
          if ((*rechit).clusterSize()>maxCLS) continue;
          //if (isMC){if (CLHEP::RandFlat::shoot(engine, 0., 100.) > chamSetEff[countTC-1]) continue;}
          testRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));        
      }
    }
    
    std::vector<int>::const_iterator it1 = idxChTraj->begin();
    std::vector<unsigned int>::const_iterator it2 = seedTypes->begin();
    std::vector<TrajectorySeed>::const_iterator it3 = seedGCM->begin();
    
    TrajectorySeed bestSeed;
    unsigned int unTypeSeed = 0;
    
    for ( ; it1 != idxChTraj->end() ; ) {
      if ( *it1 == countTC - 1 ) {
        unTypeSeed = *it2;
        bestSeed = *it3;
        break;
      }
      
      it1++;
      it2++;
      it3++;
    }
    
    if ( it1 == idxChTraj->end() ) continue;
    
    PTrajectoryStateOnDet ptsd1(bestSeed.startingState());
    DetId did(ptsd1.detId());
    const BoundPlane& bp = theService->trackingGeometry()->idToDet(did)->surface();
    TrajectoryStateOnSurface tsos = trajectoryStateTransform::transientState(ptsd1,&bp,&*theService->magneticField());
    TrajectoryStateOnSurface tsosCurrent = tsos;
    /*strKeep += TString::Format("Trajectory : (%lf, %lf, %lf) - (%lf, %lf, %lf) // \n%lf, %lf, %lf (n : %lf, %lf, %lf)\n", 
      vecSeedGP[ nIdxBest ].P1.x(), vecSeedGP[ nIdxBest ].P1.z(), vecSeedGP[ nIdxBest ].P1.y(), 
      vecSeedGP[ nIdxBest ].P2.x(), vecSeedGP[ nIdxBest ].P2.z(), vecSeedGP[ nIdxBest ].P2.y(), 
      dDirX, dDirY, dDirZ, dDirX / dDirLen, dDirY / dDirLen, dDirZ / dDirLen);*/
    //strKeep += strListRecHit;
    int nTrajHit = 0, nTrajRecHit = 0, nTestHit = 0;
    for(int c=0; c<n_ch;c++){
      GEMChamber ch = gemChambers[c];
      tsosCurrent = theService->propagator("SteppingHelixPropagatorAny")->propagate(tsosCurrent,ch.surface());
      if (!tsosCurrent.isValid()) continue;
      Global3DPoint gtrp = tsosCurrent.freeTrajectoryState()->position();
      Local3DPoint tlp = ch.surface().toLocal(gtrp);
      if (!ch.surface().bounds().inside(tlp)){continue;}
      if (ch==tch){
        int n_roll = ch.nEtaPartitions();
        double minDely = 50.;
        int mRoll = -1;
        for (int r=0;r<n_roll;r++){
          Local3DPoint rtlp = ch.etaPartition(r+1)->surface().toLocal(gtrp);
          if(minDely > abs(rtlp.y())){minDely = abs(rtlp.y()); mRoll = r+1;}
        }

        if(1 == 0 && ( mRoll == 1 || mRoll == 8 )){
          bool tester = 1;
          for (int chId = fCha; chId < lCha+1; chId++){
            if (chId == findIndex(ch.id())) continue;
            if (!firedCh[chId]) tester = 0;
          }
          if (!tester) continue;
        }

        if (mRoll == -1){cout << "no mRoll" << endl;continue;}
        int n_strip = ch.etaPartition(mRoll)->nstrips();
        double min_x = ch.etaPartition(mRoll)->centreOfStrip(0).x();
        double max_x = ch.etaPartition(mRoll)->centreOfStrip(n_strip).x();
        if ((tlp.x()>(min_x)) & (tlp.x() < (max_x)) ){
          if ( ( ( vecChamType[ countTC - 1 ] == 2 || vecChamType[ countTC - 1 ] == 1 ) && 
            ( mRoll == 1 || mRoll == 8 ) ) && unTypeSeed == 0 ) continue;
          gem_chamber_track[findIndex(ch.id())]->Fill(4.5);
          int index = findIndex(ch.id());
          double vfat = findvfat(tlp.x(), min_x, max_x);
          gem_chamber_th2D_eff[index]->Fill(vfat, mRoll);                
          gem_chamber_thxroll_eff[index]->Fill(tlp.x(), mRoll);
          gem_chamber_thxy_eff[index]->Fill(tlp.x(), gtrp.z());
          gemcrTr_g->Fill(gtrp.x(), gtrp.z(), gtrp.y());
          g_nNumTrajHit++;
          if ( nNumCurrFiredCh == 6 ) g_nNumTrajHit6++;
          
          //printf("    TrjHit : %i, (%0.5f, %0.5f, %0.5f)\n", findIndex(ch.id()), 
          //  gtrp.x(), gtrp.z(), gtrp.y());
          if ( g_dMinY > gtrp.z() ) g_dMinY = gtrp.z();
          if ( g_dMaxY < gtrp.z() ) g_dMaxY = gtrp.z();
          //printf("trajectory hit : %lf, %lf, %lf\n", gtrp.x(), gtrp.z(), gtrp.y());
          //strKeep += TString::Format("trajectory hit : %lf, %lf, %lf\n", gtrp.x(), gtrp.z(), gtrp.y());
          nTrajHit++;
          
          if ( bIsScint ) {
            gem_chamber_th2D_eff_scint[index]->Fill(vfat, mRoll);                
            gem_chamber_thxroll_eff_scint[index]->Fill(tlp.x(), mRoll);
            gem_chamber_thxy_eff_scint[index]->Fill(tlp.x(), gtrp.z());
            gemcrTrScint_g->Fill(gtrp.x(), gtrp.z(), gtrp.y());
          }
          double maxR = 99.9;
          shared_ptr<MuonTransientTrackingRecHit> tmpRecHit;
          for (auto hit : testRecHits){
            GEMDetId hitID(hit->rawId());
            if (hitID.chamberId() != ch.id()) continue;
            GlobalPoint hitGP = hit->globalPosition();
            if (abs(hitGP.x() - gtrp.x()) > maxRes) continue;
            if (abs(hitID.roll() - mRoll)>1) continue;
            double deltaR = (hitGP - gtrp).mag();
            if (maxR > deltaR){
              tmpRecHit = hit;
              maxR = deltaR;
            }
          }
          if ( !tmpRecHit && testRecHits.size() > 0 ) gem_chamber_track[findIndex(ch.id())]->Fill(6.5);
          //if(mRoll != -1)
          if(tmpRecHit){
            gem_chamber_track[findIndex(ch.id())]->Fill(5.5);
            Local3DPoint hitLP = tmpRecHit->localPosition();
            Global3DPoint recHitGP = tmpRecHit->globalPosition();
            
            gemcrCf_g->Fill(recHitGP.x(), recHitGP.z(), recHitGP.y());
          
            //printf("eff hit : %lf, %lf, %lf\n", recHitGP.x(), recHitGP.z(), recHitGP.y());
            //strKeep += TString::Format("eff hit : %lf, %lf, %lf\n", recHitGP.x(), recHitGP.z(), recHitGP.y());
            nTrajRecHit++;
            
            gem_chamber_tr2D_eff[index]->Fill(vfat, mRoll);
            gem_chamber_trxroll_eff[index]->Fill(tlp.x(), mRoll);
            gem_chamber_trxy_eff[index]->Fill(tlp.x(), gtrp.z());
            gem_chamber_local_x[index]->Fill(tlp.x(), hitLP.x());
            gem_chamber_residual[index]->Fill(tlp.x(), hitLP.x() - tlp.x());
            rh3_chamber->Fill(index);
            
            if ( bIsScint ) {
              gemcrCfScint_g->Fill(recHitGP.x(), recHitGP.z(), recHitGP.y());
              gem_chamber_tr2D_eff_scint[index]->Fill(vfat, mRoll);
              gem_chamber_trxroll_eff_scint[index]->Fill(tlp.x(), mRoll);
              gem_chamber_trxy_eff_scint[index]->Fill(tlp.x(), gtrp.z());
              gem_chamber_local_x_scint[index]->Fill(tlp.x(), hitLP.x());
              gem_chamber_residual_scint[index]->Fill(tlp.x(), hitLP.x() - tlp.x());
              rh3_chamber_scint->Fill(index);
            }
            //cout << "chamber " << index << ", Track x : " << tlp.x() <<", RecHit x : " << hitLP.x() << ", Roll : " << mRoll << endl;
            
            //printf("    recHit : %i, (%0.5f, %0.5f, %0.5f)\n", findIndex(ch.id()), 
            //  tmpRecHit->globalPosition().x(), tmpRecHit->globalPosition().z(), tmpRecHit->globalPosition().y());
            g_nNumMatched++;
          }/* else {
            nTestHit = 1;
            printf("### Yay! ###\n");
            printf(strKeep.Data());
          }*/
        }
        continue;
      }
    }
    
    //if ( nTrajRecHit != nTrajHit )
    if ( nTestHit != 0 && 1 == 0 ) {
      printf("### missing hit occurs! ###\n");
      //printf(strKeep.Data());
    }
  }
  if ( nIsLongSeed != 0 ) events_withtraj->Fill(3.5);
  if ( nIsTraceGet != 0 ) {
    if ( fChMul == 4 ) events_withtraj->Fill(4.5);
    if ( fChMul == 5 ) events_withtraj->Fill(5.5);
    if ( fChMul >= 6 ) events_withtraj->Fill(6.5);
    events_withtraj->Fill(7.5);
    
    
    //diffTrajGenRec->Fill(maxChi2, );
  }
  g_nNumTest++;
}


