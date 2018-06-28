// -*- C++ -*-
//
// Package:    L1Trigger/L1Ntuples
// Class:      L1PiXTRKTreeProducer
// 
/**\class L1PiXTRKTreeProducer L1PiXTRKTreeProducer.cc L1TriggerDPG/L1Ntuples/src/L1PiXTRKTreeProducer.cc

Description: Produce L1 Extra tree

Implementation:
     
*/
//
// Original Author:  
//         Created:  
// $Id: L1PiXTRKTreeProducer.cc,v 1.8 2012/08/29 12:44:03 jbrooke Exp $
//
//


// system include files
#include <memory>

// framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

// data formats
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

// cond formats
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"

// ROOT output stuff
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1CaloTowerDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1CaloClusterDataFormat.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "DataFormats/L1TCalorimeter/interface/CaloCluster.h"

#include "L1Trigger/L1CaloTrigger/interface/L1TkElectronTrackMatchAlgo.h"

//track trigger data formats
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
//#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
//#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

//#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"


//for pixel cluster
//
#include <map>
#include <vector>
#include <algorithm>

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Phase2ITPixelCluster/interface/Phase2ITPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
//
//for pixel cluster

//gen particles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"

#include "FastSimulation/Particle/interface/ParticleTable.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
//


// L1EG crystal
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Phase2L1CaloTrig/interface/L1EGCrystalCluster.h"
//

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// HGCAL
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerNtupleBase.h"

#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/L1THGCal/interface/HGCalCluster.h"

//
// class declaration
//

class L1PiXTRKTreeProducer : public edm::EDAnalyzer {
typedef std::vector< TTTrack < Ref_Phase2TrackerDigi_ >> L1TkTrackCollectionType;
public:
  explicit L1PiXTRKTreeProducer(const edm::ParameterSet&);
  ~L1PiXTRKTreeProducer();
  
  
private:
  virtual void beginJob(void) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

  GlobalPoint getCalorimeterPosition(double phi, double eta, double e);
  //https://github.com/truggles/L1EGRateStudies/blob/920_L1EGCrystals/src/L1EGRateStudies.cc#L143
  void doTrackMatching(const l1slhc::L1EGCrystalCluster& cluster, edm::Handle<L1TkTrackCollectionType> l1trackHandle);
  void doTrackMatching(const l1t::HGCalMulticluster& cluster, edm::Handle<L1TkTrackCollectionType> l1trackHandle);

  int getLastCopy(int eleIndx, edm::Handle<reco::GenParticleCollection>& genParticles);

  double etaToTheta(double eta);
  double thetaToEta(double theta);
  double detEtaFromEvnt(double evntEta,double z0);

  void fill(const std::vector<SimTrack>& simTracks, const std::vector<SimVertex>& simVertices );
  std::map<unsigned, unsigned> geantToIndex_;
  const double kPi = 2*asin(1.);

public:
  

private:

  // output file
  edm::Service<TFileService> fs_;

  edm::InputTag simVertex_;
  edm::InputTag simTrack_;

  edm::EDGetTokenT<vector<SimTrack> > simTrackToken_;
  edm::EDGetTokenT<vector<SimVertex> > simVertexToken_;

  
  // tree
  TTree * tree_;

  int nVtx;
  int nMeanPU;
int genpart_n;
  std::vector<float> genpart_e;
  std::vector<float> genpart_pt;
  std::vector<float> genpart_eta;
  std::vector<float> genpart_phi;
  std::vector<int> genpart_charge;
  std::vector<int> genpart_id;

  int simtrk_n;
  std::vector<float> simtrk_pt;
  std::vector<float> simtrk_eta;
  std::vector<float> simtrk_phi;
  std::vector<int> simtrk_id;
  std::vector<int> simtrk_type;
  std::vector<float> simtrk_vx;
  std::vector<float> simtrk_vy;
  std::vector<float> simtrk_vz;
  std::vector<float> sim_vx;
  std::vector<float> sim_vy;
  std::vector<float> sim_vz;

  float initialSimTkPt;
  float lastSimTkPt;
  int bremflag;
  std::vector<float> Brem_eLoss;
  std::vector<float> Brem_ptLoss;
  std::vector<float> Brempos_radius;
  std::vector<float> Brempos_x;
  std::vector<float> Brempos_y;
  std::vector<float> Brempos_z;

  std::vector<float> propgenElpart_e;
  std::vector<float> propgenElpart_pt;
  std::vector<float> propgenElpart_eta;
  std::vector<float> propgenElpart_phi;
  std::vector<int> propgenElpart_charge;
  std::vector<float> propgenElpart_x;
  std::vector<float> propgenElpart_y;
  std::vector<float> propgenElpart_z;

  std::vector<float> propgenPopart_e;
  std::vector<float> propgenPopart_pt;
  std::vector<float> propgenPopart_eta;
  std::vector<float> propgenPopart_phi;
  std::vector<int> propgenPopart_charge;
  std::vector<float> propgenPopart_x;
  std::vector<float> propgenPopart_y;
  std::vector<float> propgenPopart_z;

  int bRecHit_n;
  int fRecHit_n;
  std::vector<int> bRecHit_layer;
  std::vector<int> bRecHit_ladder;
  std::vector<int> bRecHit_module;
  std::vector<float> bRecHit_x;
  std::vector<float> bRecHit_y;
  std::vector<float> bRecHit_z;
  std::vector<float> bRecHitSize;
  std::vector<float> bRecHitSizeX;
  std::vector<float> bRecHitSizeY;
  std::vector<int> fRecHit_disk;
  std::vector<int> fRecHit_blade;
  std::vector<int> fRecHit_side;
  std::vector<int> fRecHit_panel;
  std::vector<int> fRecHit_module;
  std::vector<float> fRecHit_x;
  std::vector<float> fRecHit_y;
  std::vector<float> fRecHit_z;
  std::vector<float> fRecHitSize;
  std::vector<float> fRecHitSizeX;
  std::vector<float> fRecHitSizeY;

  int bfastsimHit_n;
  int ffastsimHit_n;
  std::vector<int> bfastsimHit_layer;
  std::vector<float> bfastsimHit_x;
  std::vector<float> bfastsimHit_y;
  std::vector<float> bfastsimHit_z;
  std::vector<int> ffastsimHit_layer;
  std::vector<float> ffastsimHit_x;
  std::vector<float> ffastsimHit_y;
  std::vector<float> ffastsimHit_z;

  int egammaC_n;
  int egammaCCluster_n;

  std::vector<float> egammaC_e;
  std::vector<float> egammaC_et;
  std::vector<float> egammaC_eta;
  std::vector<float> egammaC_phi;
  std::vector<float> egammaC_gx;
  std::vector<float> egammaC_gy;
  std::vector<float> egammaC_gz;

  std::vector<float> egammaCCluster_e;
  std::vector<float> egammaCCluster_et;
  std::vector<float> egammaCCluster_eta;
  std::vector<float> egammaCCluster_phi;
  std::vector<float> egammaCCluster_gx;
  std::vector<float> egammaCCluster_gy;
  std::vector<float> egammaCCluster_gz;
  std::vector<float> egammaCClusterP_gx;
  std::vector<float> egammaCClusterP_gy;
  std::vector<float> egammaCClusterP_gz;

  std::vector<bool> isTrackMatched;
  std::vector<float> isoConeNTrack;
  std::vector<float> isoConePtTrack;

  std::vector<float> trackHighestPt;
  std::vector<float> trackHighestPtEta;
  std::vector<float> trackHighestPtPhi;
  std::vector<float> trackHighestPtChi2;
  std::vector<float> trackHighestPtCutChi2;
  std::vector<float> trackHighestPtCutChi2Eta;
  std::vector<float> trackHighestPtCutChi2Phi;
  std::vector<float> trackHighestPtCutChi2Chi2;
  std::vector<float> trackmatchingdR;

  std::vector<bool> hgcal_isTrackMatched;
  std::vector<float> hgcal_isoConeNTrack;
  std::vector<float> hgcal_isoConePtTrack;

  std::vector<float> hgcal_trackHighestPt;
  std::vector<float> hgcal_trackHighestPtEta;
  std::vector<float> hgcal_trackHighestPtPhi;
  std::vector<float> hgcal_trackHighestPtChi2;
  std::vector<float> hgcal_trackHighestPtCutChi2;
  std::vector<float> hgcal_trackHighestPtCutChi2Eta;
  std::vector<float> hgcal_trackHighestPtCutChi2Phi;
  std::vector<float> hgcal_trackHighestPtCutChi2Chi2;
  std::vector<float> hgcal_trackmatchingdR;

  int cl3d_n_ ;
  std::vector<float> cl3d_pt_;
  std::vector<float> cl3d_energy_;
  std::vector<float> cl3d_eta_;
  std::vector<float> cl3d_phi_;
  std::vector<int> cl3d_nclu_;
  std::vector<float> cl3d_x_;
  std::vector<float> cl3d_y_;
  std::vector<int> cl3d_z_;
  //ID values
  std::vector<float> cl3d_hovere_;
  std::vector<int> cl3d_showerlength_;
  std::vector<int> cl3d_coreshowerlength_;
  std::vector<int> cl3d_firstlayer_;
  std::vector<int> cl3d_maxlayer_;
  std::vector<float> cl3d_seetot_;
  std::vector<float> cl3d_seemax_;
  std::vector<float> cl3d_spptot_;
  std::vector<float> cl3d_sppmax_;
  std::vector<float> cl3d_szz_;
  std::vector<float> cl3d_srrtot_;
  std::vector<float> cl3d_srrmax_;
  std::vector<float> cl3d_srrmean_;
  std::vector<float> cl3d_emaxe_;

  unsigned short int nEGs;
  std::vector<float> egEt;
  std::vector<float> egEta;
  std::vector<float> egPhi;
  std::vector<float> egamma_gx;
  std::vector<float> egamma_gy;
  std::vector<float> egamma_gz;

  std::vector<short int> egIEt;
  std::vector<short int> egIEta;
  std::vector<short int> egIPhi;
  std::vector<short int> egIso;
  std::vector<short int> egBx;
  std::vector<short int> egTowerIPhi;
  std::vector<short int> egTowerIEta;
  std::vector<short int> egRawEt;
  std::vector<short int> egIsoEt;
  std::vector<short int> egFootprintEt;
  std::vector<short int> egNTT;
  std::vector<short int> egShape;
  std::vector<short int> egTowerHoE;
  
  unsigned int me0SegNum;
  std::vector<uint32_t> me0SegDetId;
  std::vector<float> me0SegPosX;
  std::vector<float> me0SegPosY;
  std::vector<float> me0SegPosZ;
  std::vector<float> me0SegDirX;
  std::vector<float> me0SegDirY;
  std::vector<float> me0SegDirZ;
  std::vector<int>   me0SegNumRecHit;
  std::vector<float> me0SegDeltaPhi;

  // EDM input tags

  edm::EDGetTokenT< SiPixelRecHitCollection > tokenPixelRecHits_;
  unsigned int getLayerNumber(const DetId&, const TrackerTopology*);

  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken_;
  edm::EDGetTokenT<l1slhc::L1EGCrystalClusterCollection> egCrysClusterToken_;
  edm::EDGetTokenT<l1extra::L1EmParticleCollection> egCrysToken_;
  edm::EDGetTokenT<l1t::EGammaBxCollection> egToken_;
  edm::EDGetToken multiclusters_token_;
  edm::EDGetToken trigger_cells_token_;
  edm::EDGetToken clusters_token_;
  edm::EDGetTokenT<L1TkTrackCollectionType> L1TrackInputToken_;
  edm::EDGetTokenT<ME0SegmentCollection> me0SegmentToken_;

  edm::EDGetTokenT< std::vector< TrackingParticle > > TrackingParticleToken_;
  edm::EDGetTokenT< std::vector< TrackingVertex > > TrackingVertexToken_;

  // from HLT
  edm::ESHandle<MagneticField> magField_;
};



L1PiXTRKTreeProducer::L1PiXTRKTreeProducer(const edm::ParameterSet& iConfig):
  simVertex_ (iConfig.getParameter<edm::InputTag> ("simVertex")),
  simTrack_ (iConfig.getParameter<edm::InputTag> ("simTrack"))
{
  simVertexToken_ = consumes<vector<SimVertex> > (simVertex_);
  simTrackToken_ = consumes<vector<SimTrack> > (simTrack_);


  L1TrackInputToken_ = consumes<L1TkTrackCollectionType>(iConfig.getParameter<edm::InputTag>("L1TrackInputTag"));
  genParticleToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleToken"));
  pileupInfoToken_ = consumes<std::vector<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileupInfoToken"));
  egCrysToken_ = consumes<l1extra::L1EmParticleCollection>(iConfig.getParameter<edm::InputTag>("egCrys"));
  egCrysClusterToken_ = consumes<l1slhc::L1EGCrystalClusterCollection>(iConfig.getParameter<edm::InputTag>("egCrysCluster"));
  egToken_ = consumes<l1t::EGammaBxCollection>(iConfig.getParameter<edm::InputTag>("eg"));
  trigger_cells_token_ = consumes<l1t::HGCalTriggerCellBxCollection>(iConfig.getParameter<edm::InputTag>("TriggerCells"));
  multiclusters_token_ = consumes<l1t::HGCalMulticlusterBxCollection>(iConfig.getParameter<edm::InputTag>("Multiclusters"));
  clusters_token_ = consumes<l1t::HGCalClusterBxCollection>(iConfig.getParameter<edm::InputTag>("Clusters"));

  tokenPixelRecHits_ = consumes< SiPixelRecHitCollection >(iConfig.getParameter<edm::InputTag>("pixelRecHits"));;

  TrackingParticleToken_ = consumes< std::vector< TrackingParticle > >(iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag"));
  TrackingVertexToken_ = consumes< std::vector< TrackingVertex > >(iConfig.getParameter<edm::InputTag>("TrackingVertexInputTag"));
  
  me0SegmentToken_ = consumes<ME0SegmentCollection>(iConfig.getParameter<edm::InputTag>("me0Segment"));

  tree_=fs_->make<TTree>("L1PiXTRKTree", "L1PiXTRKTree");

  tree_->Branch("nVtx", &nVtx);
  tree_->Branch("nMeanPU", &nMeanPU);

  tree_->Branch("genPartN",          &genpart_n);
  tree_->Branch("genPartE",          &genpart_e);
  tree_->Branch("genPartPt",         &genpart_pt);
  tree_->Branch("genPartEta",        &genpart_eta);
  tree_->Branch("genPartPhi",        &genpart_phi);
  tree_->Branch("genPartCharge",     &genpart_charge);
  tree_->Branch("genPartId",         &genpart_id);
  tree_->Branch("propgenElPartE",          &propgenElpart_e);
  tree_->Branch("propgenElPartPt",         &propgenElpart_pt);
  tree_->Branch("propgenElPartEta",        &propgenElpart_eta);
  tree_->Branch("propgenElPartPhi",        &propgenElpart_phi);
  tree_->Branch("propgenElPartCharge",     &propgenElpart_charge);
  tree_->Branch("propgenElPartx",          &propgenElpart_x);
  tree_->Branch("propgenElParty",          &propgenElpart_y);
  tree_->Branch("propgenElPartz",          &propgenElpart_z);

  tree_->Branch("simTrkN", &simtrk_n);
  tree_->Branch("simTrkPt", &simtrk_pt);
  tree_->Branch("simTrkEta", &simtrk_eta);
  tree_->Branch("simTrkPhi", &simtrk_phi);
  tree_->Branch("simTrkId", &simtrk_id);
  tree_->Branch("simTrkType", &simtrk_type);
  tree_->Branch("simTrkVx", &simtrk_vx);
  tree_->Branch("simTrkVy", &simtrk_vy);
  tree_->Branch("simTrkVz", &simtrk_vz);
  tree_->Branch("simVx", &sim_vx);
  tree_->Branch("simVy", &sim_vy);
  tree_->Branch("simVz", &sim_vz);

  tree_->Branch("lastSimtkpt", &lastSimTkPt);
  tree_->Branch("initialSimtkpt", &initialSimTkPt);
  tree_->Branch("bremflag", &bremflag);
  tree_->Branch("Brempos_radius", &Brempos_radius);
  tree_->Branch("Brem_eLoss", &Brem_eLoss);
  tree_->Branch("Brem_ptLoss", &Brem_ptLoss);
  tree_->Branch("Brempos_x", &Brempos_x);
  tree_->Branch("Brempos_y", &Brempos_y);
  tree_->Branch("Brempos_z", &Brempos_z);

  tree_->Branch("propgenPoPartE",          &propgenPopart_e);
  tree_->Branch("propgenPoPartPt",         &propgenPopart_pt);
  tree_->Branch("propgenPoPartEta",        &propgenPopart_eta);
  tree_->Branch("propgenPoPartPhi",        &propgenPopart_phi);
  tree_->Branch("propgenPoPartCharge",     &propgenPopart_charge);
  tree_->Branch("propgenPoPartx",          &propgenPopart_x);
  tree_->Branch("propgenPoParty",          &propgenPopart_y);
  tree_->Branch("propgenPoPartz",          &propgenPopart_z);

  tree_->Branch("bRecHitLayer",         &bRecHit_layer);
  tree_->Branch("bRecHitLadder",        &bRecHit_ladder);
  tree_->Branch("bRecHitModule",        &bRecHit_module);
  tree_->Branch("fRecHitDisk",          &fRecHit_disk);
  tree_->Branch("fRecHitBlade",         &fRecHit_blade);
  tree_->Branch("fRecHitSide",          &fRecHit_side);
  tree_->Branch("fRecHitPanel",         &fRecHit_panel);
  tree_->Branch("fRecHitModule",        &fRecHit_module);
  tree_->Branch("bRecHitN",             &bRecHit_n);
  tree_->Branch("fRecHitN",             &fRecHit_n);
  tree_->Branch("fRecHitGx",            &fRecHit_x);
  tree_->Branch("fRecHitGy",            &fRecHit_y);
  tree_->Branch("fRecHitGz",            &fRecHit_z);
  tree_->Branch("fRhSize",           &fRecHitSize);
  tree_->Branch("fRhSizeX",          &fRecHitSizeX);
  tree_->Branch("fRhSizeY",          &fRecHitSizeY);
  tree_->Branch("bRecHitGx",            &bRecHit_x);
  tree_->Branch("bRecHitGy",            &bRecHit_y);
  tree_->Branch("bRecHitGz",            &bRecHit_z);
  tree_->Branch("bRhSize",           &bRecHitSize);
  tree_->Branch("bRhSizeX",          &bRecHitSizeX);
  tree_->Branch("bRhSizeY",          &bRecHitSizeY);

  tree_->Branch("bfastsimHitN",             &bfastsimHit_n);
  tree_->Branch("ffastsimHitN",             &ffastsimHit_n);
  tree_->Branch("bfastsimHitLayer",         &bfastsimHit_layer);
  tree_->Branch("bfastsimHitGx",            &bfastsimHit_x);
  tree_->Branch("bfastsimHitGy",            &bfastsimHit_y);
  tree_->Branch("bfastsimHitGz",            &bfastsimHit_z);
  tree_->Branch("ffastsimHitLayer",         &ffastsimHit_layer);
  tree_->Branch("ffastsimHitGx",            &ffastsimHit_x);
  tree_->Branch("ffastsimHitGy",            &ffastsimHit_y);
  tree_->Branch("ffastsimHitGz",            &ffastsimHit_z);

  tree_->Branch("egCrysN",               &egammaC_n);
  tree_->Branch("egCrysE",               &egammaC_e);
  tree_->Branch("egCrysEt",              &egammaC_et);
  tree_->Branch("egCrysEta",             &egammaC_eta);
  tree_->Branch("egCrysPhi",             &egammaC_phi);
  tree_->Branch("egCrysGx",             &egammaC_gx);
  tree_->Branch("egCrysGy",             &egammaC_gy);
  tree_->Branch("egCrysGz",             &egammaC_gz);

  tree_->Branch("egCrysClusterN",               &egammaCCluster_n);
  tree_->Branch("egCrysClusterE",               &egammaCCluster_e);
  tree_->Branch("egCrysClusterEt",              &egammaCCluster_et);
  tree_->Branch("egCrysClusterEta",             &egammaCCluster_eta);
  tree_->Branch("egCrysClusterPhi",             &egammaCCluster_phi);
  tree_->Branch("egCrysClusterGx",              &egammaCCluster_gx);
  tree_->Branch("egCrysClusterGy",              &egammaCCluster_gy);
  tree_->Branch("egCrysClusterGz",              &egammaCCluster_gz);
  tree_->Branch("egCrysClusterPGx",              &egammaCClusterP_gx);
  tree_->Branch("egCrysClusterPGy",              &egammaCClusterP_gy);
  tree_->Branch("egCrysClusterPGz",              &egammaCClusterP_gz);

  tree_->Branch("isTrackMatched", &isTrackMatched);
  tree_->Branch("isoConeNTrack", &isoConeNTrack);
  tree_->Branch("isoConePtTrack", &isoConePtTrack);

  tree_->Branch("trackHighestPt", &trackHighestPt);
  tree_->Branch("trackHighestPtEta", &trackHighestPtEta);
  tree_->Branch("trackHighestPtPhi", &trackHighestPtPhi);
  tree_->Branch("trackHighestPtChi2", &trackHighestPtChi2);
  tree_->Branch("trackHighestPtCutChi2", &trackHighestPtCutChi2);
  tree_->Branch("trackHighestPtCutChi2Eta", &trackHighestPtCutChi2Eta);
  tree_->Branch("trackHighestPtCutChi2Phi", &trackHighestPtCutChi2Phi);
  tree_->Branch("trackHighestPtCutChi2Chi2", &trackHighestPtCutChi2Chi2);
  tree_->Branch("trackmatchingdR", &trackmatchingdR);

  tree_->Branch("hgcal_isTrackMatched", &hgcal_isTrackMatched);
  tree_->Branch("hgcal_isoConeNTrack", &hgcal_isoConeNTrack);
  tree_->Branch("hgcal_isoConePtTrack", &hgcal_isoConePtTrack);

  tree_->Branch("hgcal_trackHighestPt", &hgcal_trackHighestPt);
  tree_->Branch("hgcal_trackHighestPtEta", &hgcal_trackHighestPtEta);
  tree_->Branch("hgcal_trackHighestPtPhi", &hgcal_trackHighestPtPhi);
  tree_->Branch("hgcal_trackHighestPtChi2", &hgcal_trackHighestPtChi2);
  tree_->Branch("hgcal_trackHighestPtCutChi2", &hgcal_trackHighestPtCutChi2);
  tree_->Branch("hgcal_trackHighestPtCutChi2Eta", &hgcal_trackHighestPtCutChi2Eta);
  tree_->Branch("hgcal_trackHighestPtCutChi2Phi", &hgcal_trackHighestPtCutChi2Phi);
  tree_->Branch("hgcal_trackHighestPtCutChi2Chi2", &hgcal_trackHighestPtCutChi2Chi2);
  tree_->Branch("hgcal_trackmatchingdR", &hgcal_trackmatchingdR);

  tree_->Branch("cl3d_n", &cl3d_n_, "cl3d_n/I");
  tree_->Branch("cl3d_pt", &cl3d_pt_);
  tree_->Branch("cl3d_energy", &cl3d_energy_);
  tree_->Branch("cl3d_eta", &cl3d_eta_);
  tree_->Branch("cl3d_phi", &cl3d_phi_);
  tree_->Branch("cl3d_nclu", &cl3d_nclu_);   
  tree_->Branch("cl3d_x", &cl3d_x_);
  tree_->Branch("cl3d_y", &cl3d_y_);
  tree_->Branch("cl3d_z", &cl3d_z_);   
  tree_->Branch("cl3d_hovere", &cl3d_hovere_);   
  tree_->Branch("cl3d_showerlength", &cl3d_showerlength_);
  tree_->Branch("cl3d_coreshowerlength", &cl3d_coreshowerlength_);
  tree_->Branch("cl3d_firstlayer", &cl3d_firstlayer_);
  tree_->Branch("cl3d_maxlayer", &cl3d_maxlayer_);
  tree_->Branch("cl3d_seetot", &cl3d_seetot_);
  tree_->Branch("cl3d_seemax", &cl3d_seemax_);
  tree_->Branch("cl3d_spptot", &cl3d_spptot_);
  tree_->Branch("cl3d_sppmax", &cl3d_sppmax_);
  tree_->Branch("cl3d_szz", &cl3d_szz_);
  tree_->Branch("cl3d_srrtot", &cl3d_srrtot_);
  tree_->Branch("cl3d_srrmax", &cl3d_srrmax_);
  tree_->Branch("cl3d_srrmean", &cl3d_srrmean_);
  tree_->Branch("cl3d_emaxe", &cl3d_emaxe_);

  tree_->Branch("egN",               &nEGs);
  tree_->Branch("egEt",              &egEt);
  tree_->Branch("egEta",             &egEta);
  tree_->Branch("egPhi",             &egPhi);
  tree_->Branch("egGx",             &egamma_gx);
  tree_->Branch("egGy",             &egamma_gy);
  tree_->Branch("egGz",             &egamma_gz);
  tree_->Branch("egIEt",             &egIEt);
  tree_->Branch("egIEta",             &egIEta);
  tree_->Branch("egIPhi",             &egIPhi);
  tree_->Branch("egIso",             &egIso);
  tree_->Branch("egBx",             &egBx);
  tree_->Branch("egTowerIPhi",             &egTowerIPhi);
  tree_->Branch("egTowerIEta",             &egTowerIEta);
  tree_->Branch("egRawEt",             &egRawEt);
  tree_->Branch("egIsoEt",             &egIsoEt);
  tree_->Branch("egFootprintEt",             &egFootprintEt);
  tree_->Branch("egNTT",             &egNTT);
  tree_->Branch("egShape",             &egShape);
  tree_->Branch("egTowerHoE",             &egTowerHoE);
  
  tree_->Branch("me0SegNum", &me0SegNum);
  tree_->Branch("me0SegDetId", &me0SegDetId);
  tree_->Branch("me0SegPosX", &me0SegPosX);
  tree_->Branch("me0SegPosY", &me0SegPosY);
  tree_->Branch("me0SegPosZ", &me0SegPosZ);
  tree_->Branch("me0SegDirX", &me0SegDirX);
  tree_->Branch("me0SegDirY", &me0SegDirY);
  tree_->Branch("me0SegDirZ", &me0SegDirZ);
  tree_->Branch("me0SegNumRecHit", &me0SegNumRecHit);
  tree_->Branch("me0SegDeltaPhi", &me0SegDeltaPhi);

}


L1PiXTRKTreeProducer::~L1PiXTRKTreeProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L1PiXTRKTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  genpart_n = 0;
  genpart_e.clear();
  genpart_pt.clear();
  genpart_eta.clear();
  genpart_phi.clear();
  genpart_charge.clear();
  genpart_id.clear();

  propgenElpart_e.clear();
  propgenElpart_pt.clear();
  propgenElpart_eta.clear();
  propgenElpart_phi.clear();
  propgenElpart_charge.clear();
  propgenElpart_x.clear();
  propgenElpart_y.clear();
  propgenElpart_z.clear();
  propgenPopart_e.clear();
  propgenPopart_pt.clear();
  propgenPopart_eta.clear();
  propgenPopart_phi.clear();
  propgenPopart_charge.clear();
  propgenPopart_x.clear();
  propgenPopart_y.clear();
  propgenPopart_z.clear();

  simtrk_pt.clear();
  simtrk_eta.clear();
  simtrk_phi.clear();
  simtrk_id.clear();
  simtrk_type.clear();
  simtrk_vx.clear();
  simtrk_vy.clear();
  simtrk_vz.clear();
  sim_vx.clear();
  sim_vy.clear();
  sim_vz.clear();
  
  Brempos_radius.clear();
  Brem_ptLoss.clear();
  Brem_eLoss.clear();
  Brempos_x.clear();
  Brempos_y.clear();
  Brempos_z.clear();

  bRecHit_n = 0;
  fRecHit_n = 0;
  bRecHit_layer.clear();
  bRecHit_ladder.clear();
  bRecHit_module.clear();
  fRecHit_disk.clear();
  fRecHit_blade.clear();
  fRecHit_side.clear();
  fRecHit_panel.clear();
  fRecHit_module.clear();
  bRecHit_x.clear();
  bRecHit_y.clear();
  bRecHit_z.clear();
  bRecHitSize.clear();
  bRecHitSizeX.clear();
  bRecHitSizeY.clear();
  fRecHit_x.clear();
  fRecHit_y.clear();
  fRecHit_z.clear();
  fRecHitSize.clear();
  fRecHitSizeX.clear();
  fRecHitSizeY.clear();

  bfastsimHit_n = 0;
  ffastsimHit_n = 0;
  bfastsimHit_layer.clear();
  bfastsimHit_x.clear();
  bfastsimHit_y.clear();
  bfastsimHit_z.clear();
  ffastsimHit_layer.clear();
  ffastsimHit_x.clear();
  ffastsimHit_y.clear();
  ffastsimHit_z.clear();

  egammaC_n = 0;
  egammaC_e.clear();
  egammaC_et.clear();
  egammaC_eta.clear();
  egammaC_phi.clear();
  egammaC_gx.clear();
  egammaC_gy.clear();
  egammaC_gz.clear();

  egammaCCluster_n = 0;
  egammaCCluster_e.clear();
  egammaCCluster_et.clear();
  egammaCCluster_eta.clear();
  egammaCCluster_phi.clear();
  egammaCCluster_gx.clear();
  egammaCCluster_gy.clear();
  egammaCCluster_gz.clear();
  egammaCClusterP_gx.clear();
  egammaCClusterP_gy.clear();
  egammaCClusterP_gz.clear();

  isTrackMatched.clear();
  isoConeNTrack.clear();
  isoConePtTrack.clear();

  trackHighestPt.clear();
  trackHighestPtEta.clear();
  trackHighestPtPhi.clear();
  trackHighestPtChi2.clear();
  trackHighestPtCutChi2.clear();
  trackHighestPtCutChi2Eta.clear();
  trackHighestPtCutChi2Phi.clear();
  trackHighestPtCutChi2Chi2.clear();
  trackmatchingdR.clear();

  hgcal_isTrackMatched.clear();
  hgcal_isoConeNTrack.clear();
  hgcal_isoConePtTrack.clear();

  hgcal_trackHighestPt.clear();
  hgcal_trackHighestPtEta.clear();
  hgcal_trackHighestPtPhi.clear();
  hgcal_trackHighestPtChi2.clear();
  hgcal_trackHighestPtCutChi2.clear();
  hgcal_trackHighestPtCutChi2Eta.clear();
  hgcal_trackHighestPtCutChi2Phi.clear();
  hgcal_trackHighestPtCutChi2Chi2.clear();
  hgcal_trackmatchingdR.clear();

  cl3d_n_ = 0;
  cl3d_pt_.clear();
  cl3d_energy_.clear();
  cl3d_eta_.clear();
  cl3d_phi_.clear();
  cl3d_nclu_.clear();
  cl3d_x_.clear();
  cl3d_y_.clear();
  cl3d_z_.clear();
  cl3d_hovere_.clear();
  cl3d_showerlength_.clear();
  cl3d_coreshowerlength_.clear();
  cl3d_firstlayer_.clear();
  cl3d_maxlayer_.clear();
  cl3d_seetot_.clear();
  cl3d_seemax_.clear();
  cl3d_spptot_.clear();
  cl3d_sppmax_.clear();
  cl3d_szz_.clear();
  cl3d_srrtot_.clear();
  cl3d_srrmax_.clear();
  cl3d_srrmean_.clear();
  cl3d_emaxe_.clear();

  nEGs = 0;
  egEt.clear();
  egEta.clear();
  egPhi.clear();
  egamma_gx.clear();
  egamma_gy.clear();
  egamma_gz.clear();
  egIEt.clear();
  egIEta.clear();
  egIPhi.clear();
  egIso.clear();
  egBx.clear();
  egTowerIPhi.clear();
  egTowerIEta.clear();
  egRawEt.clear();
  egIsoEt.clear();
  egFootprintEt.clear();
  egNTT.clear();
  egShape.clear();
  egTowerHoE.clear();
  
  me0SegNum = 0;
  me0SegDetId.clear();
  me0SegPosX.clear();
  me0SegPosY.clear();
  me0SegPosZ.clear();
  me0SegDirX.clear();
  me0SegDirY.clear();
  me0SegDirZ.clear();
  me0SegNumRecHit.clear();
  me0SegDeltaPhi.clear();

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);

  edm::Handle<SiPixelRecHitCollection> recHits;
  iEvent.getByToken( tokenPixelRecHits_, recHits );

  edm::ESHandle< TrackerGeometry > geomHandle;
  iSetup.get< TrackerDigiGeometryRecord >().get(geomHandle);
  const TrackerGeometry* tkGeom = &(*geomHandle);

  edm::ESHandle< TrackerTopology > tTopoHandle;
  iSetup.get< TrackerTopologyRcd >().get(tTopoHandle);
  const TrackerTopology* tTopo = tTopoHandle.product();

  edm::Handle<std::vector<PileupSummaryInfo>> puInfoCollection;
  iEvent.getByToken(pileupInfoToken_, puInfoCollection);

  if (!puInfoCollection.isValid()) {
    throw cms::Exception("ProductNotValid") << "pileupInfoSource not valid";
  }

  // Loop over vector, find in-time entry, then store the relevant info
  std::vector<PileupSummaryInfo>::const_iterator puItr = puInfoCollection->begin();
  std::vector<PileupSummaryInfo>::const_iterator puEnd = puInfoCollection->end();
  for( ; puItr != puEnd; ++puItr) {
    int bx = puItr->getBunchCrossing();
    if (bx == 0) {
      nMeanPU = puItr->getTrueNumInteractions();
      nVtx    = puItr->getPU_NumInteractions();
      break;
    }
  }


  genpart_n = genParticles->size();
  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const reco::GenParticle & p = (*genParticles)[i];

    genpart_e.push_back(      p.energy() );
    genpart_pt.push_back(     p.pt()     );
    genpart_eta.push_back(    p.eta()    );
    genpart_phi.push_back(    p.phi()    );
    genpart_charge.push_back( p.charge() );
    genpart_id.push_back(     p.pdgId()  );
  }

  //cout << "find gen electron/ positron" << endl;
  //// find the last electron 
  //int eleIndx = -1;
  //std::vector<int> pids(1,11);

  //for(size_t i = 0; i < genParticles->size(); ++ i) {
  //  const reco::GenParticle & part = (*genParticles)[i];
  //  if( part.numberOfMothers() != 1 ) continue;
  //  if( (part.status()>=20 && part.status()<=29) || part.mother(0)->pdgId() == 23 ){ //

  //    if(std::find(pids.begin(),pids.end(),part.pdgId())!=pids.end()){
  //      eleIndx = i;
  //    }
  //  }
  //}

  //cout << "ele index bf getLastCopy: " << eleIndx << endl;
  //int eleIndx_last = getLastCopy(eleIndx, genParticles);
 
  //cout << "ele index: " << eleIndx_last << endl;

  //// find the last positron 
  //int posIndx = -1;
  //std::vector<int> pids_(1,-11);

  //for(size_t i = 0; i < genParticles->size(); ++ i) {
  //  const reco::GenParticle & part = (*genParticles)[i];
  //  if( part.numberOfMothers() != 1 ) continue;
  //  if( (part.status()>=20 && part.status()<=29) || part.mother(0)->pdgId() == 23 ){ //

  //    if(std::find(pids_.begin(),pids_.end(),part.pdgId())!=pids_.end()){
  //      posIndx = i;
  //    }
  //  }
  //}

  //int posIndx_last = getLastCopy(posIndx, genParticles);
  //cout << "pos index: " << posIndx_last << endl;

  //std::vector<int> dilep_index;
  //dilep_index.push_back(eleIndx_last);
  //dilep_index.push_back(posIndx_last);

  //for(size_t i = 0; i < dilep_index.size(); ++ i) {
  //  const reco::GenParticle & p = (*genParticles)[dilep_index[i]];

  //  genpart_e.push_back(      p.energy() );
  //  genpart_pt.push_back(     p.pt()     );
  //  genpart_eta.push_back(    p.eta()    );
  //  genpart_phi.push_back(    p.phi()    );
  //  genpart_charge.push_back( p.charge() );
  //  genpart_id.push_back(     p.pdgId()  );
  //}

  // simtrack info
  edm::Handle< edm::SimTrackContainer >   simTrackHandle;
  edm::Handle< edm::SimVertexContainer >  simVertex;
  iEvent.getByToken( simVertexToken_, simVertex);
  iEvent.getByToken( simTrackToken_, simTrackHandle );

  edm::SimTrackContainer::const_iterator iterSimTracks;
  simtrk_n = 0;
  for (iterSimTracks = simTrackHandle->begin(); iterSimTracks != simTrackHandle->end(); ++iterSimTracks) {
    simtrk_n++;
    simtrk_pt.push_back( iterSimTracks->momentum().pt() );
    simtrk_eta.push_back( iterSimTracks->momentum().eta() );
    simtrk_phi.push_back( iterSimTracks->momentum().phi() );
    simtrk_id.push_back( iterSimTracks->trackId() );
    simtrk_type.push_back( iterSimTracks->type() );
    int index = iterSimTracks->vertIndex();
    simtrk_vx.push_back( simVertex->at(index).position().x() );
    simtrk_vy.push_back( simVertex->at(index).position().y() );
    simtrk_vz.push_back( simVertex->at(index).position().z() );
    //std::cout << simtrk_n << " th track, " << " vertex process id: " << simVertex->at(index).processType() << " track id: " << iterSimTracks->trackId()<< std::endl;
    //cout << "simtrack pt: " << iterSimTracks->momentum().pt() << " energy: " << iterSimTracks->momentum().E() << ", eta: " << iterSimTracks->momentum().eta() << ", phi: " << iterSimTracks->momentum().phi() << ", type: " << iterSimTracks->type() << ", simvertex: " << simVertex->at(index).position().z() << " parent track id: " << simVertex->at(index). parentIndex() << endl;
    // index==0 gets the primary vertex;
    if(index==0){
      sim_vx.push_back( simVertex->at(index).position().x() );
      sim_vy.push_back( simVertex->at(index).position().y() );
      sim_vz.push_back( simVertex->at(index).position().z() );
    }
  }


  // brem information 
  //if(simTrackHandle.isValid())
  if(simtrk_n != 0){
   // Find Brem vertex with SimTrack & SimVertex ( copy & paste RecoEgamma/EgammaMCTools/src/ElectronMCTruthFinder.cc 
   edm::SimTrackContainer theSimTracks = *simTrackHandle.product(); // in /SimDataFormats/Track/interface/SimTrackContainer.h, typedef std::vector<SimTrack> SimTrackContainer;
   edm::SimVertexContainer theSimVertices = *simVertex.product();

   std::vector<SimTrack> electronTracks;
   SimVertex primVtx;

   fill(theSimTracks,  theSimVertices);

  int iPV=-1;
  //int partType1=0;
  //int partType2=0;
  std::vector<SimTrack>::const_iterator iFirstSimTk = theSimTracks.begin();
  if (  !(*iFirstSimTk).noVertex() ) {
    iPV =  (*iFirstSimTk).vertIndex();

    int vtxId =   (*iFirstSimTk).vertIndex();
    primVtx = theSimVertices[vtxId];

    //partType1 = (*iFirstSimTk).type();
  } else {
    //std::cout << " First track has no vertex " << std::endl;
  }

  // collect electron tracks
  int npv=0;
  for (std::vector<SimTrack>::const_iterator iSimTk = theSimTracks.begin(); iSimTk != theSimTracks.end(); ++iSimTk){
    if (  (*iSimTk).noVertex() ) continue;

    //int vertexId = (*iSimTk).vertIndex();
    //SimVertex vertex = theSimVertices[vertexId];

    //std::cout << " Particle type " <<  (*iSimTk).type() << " Sim Track ID " << (*iSimTk).trackId() << " momentum " << (*iSimTk).momentum() <<  " vertex position " << vertex.position() << " vertex ID " << vertexId  << std::endl;  
    if ( (*iSimTk).vertIndex() == iPV ) {
      npv++;
      if ( (*iSimTk).type() == 11) { // consider only electron track
        //std::cout << " Found a primary electron with ID  " << (*iSimTk).trackId() << " momentum " << (*iSimTk).momentum() <<  std::endl;
        electronTracks.push_back( *iSimTk );
      }
    }
  }

  bremflag = -1;
  if(electronTracks.size() != 0){
      bremflag = 0;
      //cout << "electron track found... energy: " << electronTracks[0].momentum().e() << endl;
   
      //SimTrack trLast = theSimTracks[0]; 
      SimTrack trLast = electronTracks[0]; 
      //unsigned int eleId = theSimTracks[0].trackId();
      unsigned int eleId = electronTracks[0].trackId();
      float remainingPt =trLast.momentum().pt();
      float remainingEnergy =trLast.momentum().e();

      //math::XYZTLorentzVectorD motherMomentum((theSimTracks[0]).momentum().x(),
      //                                        (theSimTracks[0]).momentum().y(),
      //                                        (theSimTracks[0]).momentum().z(),
      //                                        (theSimTracks[0]).momentum().e());

      //int eleVtxIndex= (theSimTracks[0]).vertIndex();

      //int iPV =  (theSimTracks[0]).vertIndex();

      initialSimTkPt = (electronTracks[0]).momentum().pt();
      lastSimTkPt = 0.;
      for (std::vector<SimTrack>::const_iterator iSimTk = theSimTracks.begin(); iSimTk != theSimTracks.end(); ++iSimTk){

        if (  (*iSimTk).noVertex() )                    continue;
        if ( (*iSimTk).vertIndex() == iPV )             continue;

        int vertexId1 = (*iSimTk).vertIndex();
        SimVertex vertex1 = theSimVertices[vertexId1];
        int vertexId2 = trLast.vertIndex();
        int motherId=-1;

        if(  (  vertexId1 ==  vertexId2 ) && ( (*iSimTk).type() == (electronTracks[0]).type() ) && trLast.type() == 22   ) {

          float ptLoss = (remainingPt - ( (*iSimTk).momentum() + trLast.momentum()).pt()) / electronTracks[0].momentum().pt();
          float eLoss = (remainingEnergy - ( (*iSimTk).momentum() + trLast.momentum()).e()) / electronTracks[0].momentum().e();

          //std::cout << " Found a brem vertex track ID  " << (*iSimTk).trackId() << " momentum pt " << (*iSimTk).momentum().pt() << " energy " << (*iSimTk).momentum().e() << std::endl;
          //std::cout << " photon track ID  " << trLast.trackId() << " momentum pt " << trLast.momentum().pt() << " energy " << trLast.momentum().e() << std::endl;
          //std::cout << "ptLoss: " << ptLoss << " eLoss: " << eLoss << std::endl;

          if ( vertex1.parentIndex()  ) {

            unsigned  motherGeantId = vertex1.parentIndex();
            std::map<unsigned, unsigned >::iterator association = geantToIndex_.find( motherGeantId );
            if(association != geantToIndex_.end() )
              motherId = association->second;

            if ( theSimTracks[motherId].trackId() == eleId ) {
              //std::cout << "process type: " << simVertex->at(vertexId1).processType() << std::endl;
              eleId= (*iSimTk).trackId();
              remainingPt = (*iSimTk).momentum().pt();
              remainingEnergy = (*iSimTk).momentum().e();

              lastSimTkPt = (*iSimTk).momentum().pt();
              bremflag = 1;
              Brempos_radius.push_back(sqrt(pow(simVertex->at(vertexId1).position().x(),2)+pow(simVertex->at(vertexId1).position().y(),2)));
              Brempos_x.push_back(simVertex->at(vertexId1).position().x());
              Brempos_y.push_back(simVertex->at(vertexId1).position().y());
              Brempos_z.push_back(simVertex->at(vertexId1).position().z());
              Brem_eLoss.push_back(eLoss);
              Brem_ptLoss.push_back(ptLoss);

              //cout << " brem x: " << simVertex->at(vertexId1).position().x() << " y: " << simVertex->at(vertexId1).position().y() << " z: " << simVertex->at(vertexId1).position().z() << endl;
            }
          } else {
            //std::cout << " This vertex has no parent tracks " <<  std::endl;
          }
        }
        trLast=(*iSimTk);
      }
      
      //
      if(!bremflag){
       //lastSimTkPt = (theSimTracks[0]).momentum().pt();
       lastSimTkPt = (electronTracks[0]).momentum().pt();
      }
   }
   else{
      cout << "electron track not found... " << endl;
      //const reco::GenParticle & p = (*genParticles)[dilep_index[0]];
      const reco::GenParticle & p = (*genParticles)[0];

      cout << "gen electron pt: " << p.pt() << " eta: " << p.eta() << " phi: " << p.phi() << endl;
   }
  }

  //cout << "brem info stored..." << endl;

  // tracking particles
  edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  edm::Handle< std::vector< TrackingVertex > > TrackingVertexHandle;
  iEvent.getByToken(TrackingParticleToken_, TrackingParticleHandle);
  iEvent.getByToken(TrackingVertexToken_, TrackingVertexHandle);

  // check the indix of gen particle; 
  // Propagated true electron to ECAL
  reco::GenParticleCollection genParticle = *genParticles.product();
  reco::Candidate::PolarLorentzVector trueElectron;
  RawParticle electron(genParticle[0].p4());
  //RawParticle electron(genParticle[eleIndx_last].p4());
  electron.setVertex(genParticle[0].vertex().x(), genParticle[0].vertex().y(), genParticle[0].vertex().z(), 0.);
  //electron.setVertex(genParticle[eleIndx_last].vertex().x(), genParticle[eleIndx_last].vertex().y(), genParticle[eleIndx_last].vertex().z(), 0.);
  electron.setCharge(-1.);
  //electron.setCharge(genParticle[eleIndx_last].charge());
  electron.setMass(0.000511);

  BaseParticlePropagator prop_electron(electron, 0., 0., 4.); // 4. is the magnetic field strength (not sure why they expect hardcoded homogenous field)
  BaseParticlePropagator start_el(prop_electron);
  prop_electron.propagateToEcalEntrance();

  if(prop_electron.getSuccess()!=0)
  {
    trueElectron = reco::Candidate::PolarLorentzVector(prop_electron.E()*sin(prop_electron.vertex().theta()), prop_electron.vertex().eta(), prop_electron.vertex().phi(), 0.);
  }
  else
  {
    // something failed
    trueElectron = genParticle[0].polarP4();
    //trueElectron = genParticle[eleIndx_last].polarP4();
  }

  propgenElpart_e.push_back(trueElectron.energy());
  propgenElpart_pt.push_back(     trueElectron.pt()     );
  propgenElpart_eta.push_back(    trueElectron.eta()    );
  propgenElpart_phi.push_back(    trueElectron.phi()    );
  propgenElpart_charge.push_back( electron.charge() );
  GlobalPoint lcpos= getCalorimeterPosition(trueElectron.phi(), trueElectron.eta(), trueElectron.energy());
  propgenElpart_x.push_back(lcpos.x());
  propgenElpart_y.push_back(lcpos.y());
  propgenElpart_z.push_back(lcpos.z());


  // test for fast simulation for an electrons propagate through BPIX1
  //reco::Candidate::PolarLorentzVector testElectron;

  double PXB_r[4] = { 2.90, 7.0146, 11.7753, 15.7388};
  double PXB_z_max = 19.9700;

  //std::cout << "==============================Gen electron information============================================" << std::endl;
  //std::cout << "pt: " << genParticle[0].pt() << " eta: " << genParticle[0].eta() << std::endl;
  //std::cout << "initial r: " << sqrt(pow(genParticle[0].vertex().x(),2.) + pow(genParticle[0].vertex().y(),2.)) << " z: " << genParticle[0].vertex().z() << std::endl;

  for( int i = 0; i < 4; i++){
     RawParticle electron_test(genParticle[0].p4());
     electron_test.setVertex(genParticle[0].vertex().x(), genParticle[0].vertex().y(), genParticle[0].vertex().z(), 0.);
     electron_test.setCharge(-1.);
     electron_test.setMass(0.000511);

     BaseParticlePropagator prop_electron_test(electron_test, 0., 0., 4.); 
     prop_electron_test.setPropagationConditions(PXB_r[i], PXB_z_max, true); // PXB geometry from http://cms-tklayout.web.cern.ch/cms-tklayout/layouts/cmssw-models/I_OT613_200_IT4025/layoutpixel.html, if the last term is set true, and if the particle does not cross the barrel retrun 0
     bool done = prop_electron_test.propagate(); // vertex will be updated to the location of intersection
     
     if(done){
       if(fabs(prop_electron_test.getSuccess()) == 1){ // intersection is in the barrel
         bfastsimHit_n++;
         bfastsimHit_layer.push_back(i + 1);
         bfastsimHit_x.push_back(prop_electron_test.vertex().x());
         bfastsimHit_y.push_back(prop_electron_test.vertex().y());
         bfastsimHit_z.push_back(prop_electron_test.vertex().z());

         //std::cout << "Intersection point information at PXB " << i + 1 << std::endl;
         //testElectron = reco::Candidate::PolarLorentzVector(prop_electron_test.E()*sin(prop_electron_test.vertex().theta()), prop_electron_test.vertex().eta(), prop_electron_test.vertex().phi(), 0.); 
         //std::cout << "propagation done succesfully, getSuccess: " << prop_electron_test.getSuccess() << " r: " << sqrt(pow(prop_electron_test.vertex().x(),2.)+pow(prop_electron_test.vertex().y(), 2.)) << " z: " << prop_electron_test.vertex().z()  << std::endl;
       }
     }
  }

  double FPIX_1_r_max = 16.0839;
  double FPIX_1_r_min = 2.90;
  double FPIX_1_z[8] = {25.0000, 31.9760, 40.8990, 52.3110, 66.9080, 85.5780, 109.4570, 140.0000};

  for( int i = 0; i < 8; i++){
     RawParticle electron_test(genParticle[0].p4());
     electron_test.setVertex(genParticle[0].vertex().x(), genParticle[0].vertex().y(), genParticle[0].vertex().z(), 0.);
     electron_test.setCharge(-1.);
     electron_test.setMass(0.000511);

     BaseParticlePropagator prop_electron_test(electron_test, 0., 0., 4.);
     prop_electron_test.setPropagationConditions(FPIX_1_r_max, FPIX_1_z[i], true); // PXB geometry from http://cms-tklayout.web.cern.ch/cms-tklayout/layouts/cmssw-models/I_OT613_200_IT4025/layoutpixel.html, if the last term is set true, and if the particle does not cross the barrel retrun 0
     bool done = prop_electron_test.propagate(); // vertex will be updated to the location of intersection

     if(done){
       if(fabs(prop_electron_test.getSuccess()) == 2){
         ffastsimHit_n++;
         ffastsimHit_layer.push_back(i + 1);
         ffastsimHit_x.push_back(prop_electron_test.vertex().x());
         ffastsimHit_y.push_back(prop_electron_test.vertex().y());
         ffastsimHit_z.push_back(prop_electron_test.vertex().z());

         //std::cout << "Intersection point information at FPIX_1 " << i + 1 << std::endl;
         //testElectron = reco::Candidate::PolarLorentzVector(prop_electron_test.E()*sin(prop_electron_test.vertex().theta()), prop_electron_test.vertex().eta(), prop_electron_test.vertex().phi(), 0.); 
         //std::cout << "propagation done succesfully, getSuccess: " << prop_electron_test.getSuccess() << " r: " << sqrt(pow(prop_electron_test.vertex().x(),2.)+pow(prop_electron_test.vertex().y(), 2.)) << " z: " << prop_electron_test.vertex().z()  << std::endl;
       }
     }
  }
  
  // Putting the infos for ME0
  edm::Handle<ME0SegmentCollection> me0SegmentCollection;
  iEvent.getByToken(me0SegmentToken_, me0SegmentCollection);
  
  edm::ESHandle<ME0Geometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  const ME0Geometry* ME0Geometry_ = &*hGeom;

  for ( ME0SegmentCollection::const_iterator it = me0SegmentCollection->begin() ; 
    it != me0SegmentCollection->end() ; it++ )
  {
    ME0DetId me0Id = it->me0DetId();
    LocalPoint lpCurr = it->localPosition();
    LocalVector ldCurr = it->localDirection();
    
    GlobalPoint gpCurr = ME0Geometry_->idToDet(me0Id)->surface().toGlobal(lpCurr);
    GlobalVector gdCurr = ME0Geometry_->idToDet(me0Id)->surface().toGlobal(lpCurr + ldCurr) - gpCurr;
    
    me0SegDetId.push_back(me0Id.rawId());
    
    me0SegPosX.push_back(gpCurr.x());
    me0SegPosY.push_back(gpCurr.y());
    me0SegPosZ.push_back(gpCurr.z());
    me0SegDirX.push_back(gdCurr.x());
    me0SegDirY.push_back(gdCurr.y());
    me0SegDirZ.push_back(gdCurr.z());
    
    me0SegNumRecHit.push_back(it->nRecHits());
    me0SegDeltaPhi.push_back(it->deltaPhi());
    
    printf("(%lf, %lf, %lf) - (%lf, %lf, %lf); %i, %lf\n", 
      me0SegPosX[ me0SegNum ], me0SegPosY[ me0SegNum ], me0SegPosZ[ me0SegNum ], 
      me0SegDirX[ me0SegNum ], me0SegDirY[ me0SegNum ], me0SegDirZ[ me0SegNum ], 
      me0SegNumRecHit[ me0SegNum ], me0SegDeltaPhi[ me0SegNum ]);
    
    me0SegNum++;
  }
  printf("# of ME0 Segments : %i\n", me0SegNum);


  // Propagated true positron to ECAL
  //reco::Candidate::PolarLorentzVector truePositron;

  //RawParticle positron(genParticle[posIndx_last].p4());
  //positron.setVertex(genParticle[posIndx_last].vertex().x(), genParticle[posIndx_last].vertex().y(), genParticle[posIndx_last].vertex().z(), 0.);
  //positron.setCharge(genParticle[posIndx_last].charge());
  //positron.setMass(0.000511);

  //BaseParticlePropagator prop_positron(positron, 0., 0., 4.); // 4. is the magnetic field strength (not sure why they expect hardcoded homogenous field)
  //BaseParticlePropagator start_po(prop_positron);
  //prop_positron.propagateToEcalEntrance();

  //if(prop_positron.getSuccess()!=0)
  //{
  //  truePositron = reco::Candidate::PolarLorentzVector(prop_positron.E()*sin(prop_positron.vertex().theta()), prop_positron.vertex().eta(), prop_positron.vertex().phi(), 0.);
  //}
  //else
  //{
  //  // something failed
  //  //truePositron = genParticle[0].polarP4();
  //  truePositron = genParticle[posIndx_last].polarP4();
  //}

  //propgenPopart_e.push_back(truePositron.energy());
  //propgenPopart_pt.push_back(     truePositron.pt()     );
  //propgenPopart_eta.push_back(    truePositron.eta()    );
  //propgenPopart_phi.push_back(    truePositron.phi()    );
  //propgenPopart_charge.push_back( positron.charge() );
  //lcpos= getCalorimeterPosition(truePositron.phi(), truePositron.eta(), truePositron.energy());
  //propgenPopart_x.push_back(lcpos.x());
  //propgenPopart_y.push_back(lcpos.y());
  //propgenPopart_z.push_back(lcpos.z());

  //cout << "propagated gen info stored..." << endl;

  SiPixelRecHitCollection::const_iterator detUnitIt    = recHits->begin();
  SiPixelRecHitCollection::const_iterator detUnitItEnd = recHits->end();
  for ( ; detUnitIt != detUnitItEnd; detUnitIt++ ) {
      DetId detId = DetId(detUnitIt->detId());
      int subid = detId.subdetId();
      SiPixelRecHitCollection::DetSet::const_iterator recHitIt    = detUnitIt->begin();
      SiPixelRecHitCollection::DetSet::const_iterator recHitItEnd = detUnitIt->end();
      for ( ; recHitIt != recHitItEnd; ++recHitIt) {
          LocalPoint  lp = recHitIt->localPosition();
          GlobalPoint gp = ( (geomHandle.product())->idToDet(detId) )->surface().toGlobal(lp);
          SiPixelRecHit::ClusterRef const& Cluster = recHitIt->cluster();
          if ( (gp.perp() < 20 && fabs(gp.z()) < 150) || (gp.perp() < 25 && fabs(gp.z()) > 150) ) { // drop outer tracker 
                if ( subid==PixelSubdetector::PixelEndcap ){
                     //std::cout << "rechit endcap" << " x: " << gp.x() << " y: " << gp.y() << " z: " << gp.z() << std::endl;
                     fRecHit_n ++;
                     fRecHit_disk.push_back(  tTopo->pxfDisk(detId)    );
                     fRecHit_blade.push_back( tTopo->pxfBlade(detId)    );
                     fRecHit_side.push_back(  tTopo->pxfSide(detId)    );
                     fRecHit_panel.push_back(  tTopo->pxfPanel(detId)    );
                     fRecHit_module.push_back( tTopo->pxfModule(detId)    );
                     fRecHit_x.push_back(    gp.x() );
                     fRecHit_y.push_back(    gp.y() );
                     fRecHit_z.push_back(    gp.z() );
                     fRecHitSize.push_back(   Cluster->size()  );
                     fRecHitSizeX.push_back(  Cluster->sizeX() );
                     fRecHitSizeY.push_back(  Cluster->sizeY() );
                }
                if ( subid==PixelSubdetector::PixelBarrel ){
                     //std::cout << "rechit Barrel" << " x: " << gp.x() << " y: " << gp.y() << " z: " << gp.z() << std::endl;
                     bRecHit_n ++;
                     bRecHit_layer.push_back(  getLayerNumber(detId,tTopo ) );
                     bRecHit_ladder.push_back( tTopo->pxbLadder(detId) );
                     bRecHit_module.push_back( tTopo->pxbModule(detId) );
                     bRecHit_x.push_back(    gp.x() );
                     bRecHit_y.push_back(    gp.y() );
                     bRecHit_z.push_back(    gp.z() );
                     bRecHitSize.push_back(   Cluster->size()  );
                     bRecHitSizeX.push_back(  Cluster->sizeX() );
                     bRecHitSizeY.push_back(  Cluster->sizeY() );
                }
          }
      } // close recHits loop
  } // close detUnits loop

  edm::Handle< l1extra::L1EmParticleCollection > EgammaC;
  iEvent.getByToken( egCrysToken_, EgammaC );
  //egammaC_n = EgammaC->size();
  for(l1extra::L1EmParticleCollection::const_iterator it = EgammaC->begin(); it!=EgammaC->end(); ++it){
    //std::cout << "sc pt: " << it->et() << " eta: " << it->eta() << " phi: " << it->phi() << std::endl;
    egammaC_n ++;
    egammaC_e.push_back(it->energy());
    egammaC_et.push_back(it->et());
    egammaC_eta.push_back(it->eta());
    egammaC_phi.push_back(it->phi());
    GlobalPoint pos= getCalorimeterPosition(it->phi(), it->eta(), it->energy());
    egammaC_gx.push_back(pos.x());
    egammaC_gy.push_back(pos.y());
    egammaC_gz.push_back(pos.z());
  }

  edm::Handle< l1slhc::L1EGCrystalClusterCollection > EgammaCCluster;
  iEvent.getByToken( egCrysClusterToken_, EgammaCCluster );

  edm::Handle<L1TkTrackCollectionType> l1trackHandle;
  iEvent.getByToken(L1TrackInputToken_, l1trackHandle);

  for(l1slhc::L1EGCrystalClusterCollection::const_iterator it = EgammaCCluster->begin(); it!=EgammaCCluster->end(); ++it){
    //std::cout << "sc pt: " << it->et() << " eta: " << it->eta() << " phi: " << it->phi() << std::endl;
    egammaCCluster_n ++;
    egammaCCluster_e.push_back(it->energy());
    egammaCCluster_et.push_back(it->et());
    egammaCCluster_eta.push_back(it->eta());
    egammaCCluster_phi.push_back(it->phi());
    GlobalPoint pos= getCalorimeterPosition(it->phi(), it->eta(), it->energy());
    egammaCCluster_gx.push_back(pos.x());
    egammaCCluster_gy.push_back(pos.y());
    egammaCCluster_gz.push_back(pos.z());
    GlobalVector position = it->position();;
    egammaCClusterP_gx.push_back(position.x());
    egammaCClusterP_gy.push_back(position.y());
    egammaCClusterP_gz.push_back(position.z());

    doTrackMatching(*it, l1trackHandle);
  }
 

  // retrieve trigger cells
  edm::Handle<l1t::HGCalTriggerCellBxCollection> trigger_cells_h;
  iEvent.getByToken(trigger_cells_token_, trigger_cells_h);
  const l1t::HGCalTriggerCellBxCollection& trigger_cells = *trigger_cells_h;


  edm::Handle<l1t::HGCalClusterBxCollection> clusters_h;
  iEvent.getByToken(clusters_token_, clusters_h);
  const l1t::HGCalClusterBxCollection& clusters = *clusters_h;


  // HGCAL 3D cluster
  edm::Handle<l1t::HGCalMulticlusterBxCollection> multiclusters_h;
  iEvent.getByToken(multiclusters_token_, multiclusters_h);
  const l1t::HGCalMulticlusterBxCollection& multiclusters = *multiclusters_h;

  // retrieve geometry
  edm::ESHandle<HGCalTriggerGeometryBase> geometry;
  iSetup.get<IdealGeometryRecord>().get(geometry);

  for(auto tc_itr=trigger_cells.begin(0); tc_itr!=trigger_cells.end(0); tc_itr++)
  {
    if(tc_itr->hwPt()>0)
    {
      //std::cout << "tc z: " << tc_itr->position().z() << std::endl;
    }
  }

  for(auto cl_itr=clusters.begin(0); cl_itr!=clusters.end(0); cl_itr++)
  {
    //std::cout << "2D cluster z: " << cl_itr->position().z() << std::endl;
  }

  for(auto cl3d_itr=multiclusters.begin(0); cl3d_itr!=multiclusters.end(0); cl3d_itr++)
  {
    cl3d_n_++;
    // physical values 
    cl3d_pt_.emplace_back(cl3d_itr->pt());
    cl3d_energy_.emplace_back(cl3d_itr->energy());
    cl3d_eta_.emplace_back(cl3d_itr->eta());
    cl3d_phi_.emplace_back(cl3d_itr->phi());
    cl3d_nclu_.emplace_back(cl3d_itr->constituents().size());


    cl3d_hovere_.emplace_back(cl3d_itr->hOverE());
    cl3d_showerlength_.emplace_back(cl3d_itr->showerLength());
    cl3d_coreshowerlength_.emplace_back(cl3d_itr->coreShowerLength());
    cl3d_firstlayer_.emplace_back(cl3d_itr->firstLayer());
    cl3d_maxlayer_.emplace_back(cl3d_itr->maxLayer());
    cl3d_seetot_.emplace_back(cl3d_itr->sigmaEtaEtaTot());
    cl3d_seemax_.emplace_back(cl3d_itr->sigmaEtaEtaMax());
    cl3d_spptot_.emplace_back(cl3d_itr->sigmaPhiPhiTot());
    cl3d_sppmax_.emplace_back(cl3d_itr->sigmaPhiPhiMax());
    cl3d_szz_.emplace_back(cl3d_itr->sigmaZZ());
    cl3d_srrtot_.emplace_back(cl3d_itr->sigmaRRTot());
    cl3d_srrmax_.emplace_back(cl3d_itr->sigmaRRMax());
    cl3d_srrmean_.emplace_back(cl3d_itr->sigmaRRMean());
    cl3d_emaxe_.emplace_back(cl3d_itr->eMax()/cl3d_itr->energy());


    GlobalPoint xyz = cl3d_itr->position();
    cl3d_x_.emplace_back(xyz.x());
    cl3d_y_.emplace_back(xyz.y());
    cl3d_z_.emplace_back(xyz.z());

    doTrackMatching(*cl3d_itr, l1trackHandle);
  }


  edm::Handle<l1t::EGammaBxCollection> eg;
  iEvent.getByToken(egToken_,   eg);
  unsigned maxL1Upgrade = 60;

  for (int ibx = eg->getFirstBX(); ibx <= eg->getLastBX(); ++ibx) {
    for (l1t::EGammaBxCollection::const_iterator it=eg->begin(ibx); it!=eg->end(ibx) && nEGs<maxL1Upgrade; it++){
      if (it->pt() > 0){
        egEt .push_back(it->pt());
        egEta.push_back(it->eta());
        egPhi.push_back(it->phi());
        GlobalPoint pos= getCalorimeterPosition(it->phi(), it->eta(), it->energy());
        egamma_gx.push_back(pos.x());
        egamma_gy.push_back(pos.y());
        egamma_gz.push_back(pos.z());
        egIEt .push_back(it->hwPt());
        egIEta.push_back(it->hwEta());
        egIPhi.push_back(it->hwPhi());
        egIso.push_back(it->hwIso());
        egBx .push_back(ibx);
        egTowerIPhi.push_back(it->towerIPhi());
        egTowerIEta.push_back(it->towerIEta());
        egRawEt.push_back(it->rawEt());
        egIsoEt.push_back(it->isoEt());
        egFootprintEt.push_back(it->footprintEt());
        egNTT.push_back(it->nTT());
        egShape.push_back(it->shape());
        egTowerHoE.push_back(it->towerHoE());
        nEGs++;
      }
    }
  }

  tree_->Fill();

}

GlobalPoint L1PiXTRKTreeProducer::getCalorimeterPosition(double phi, double eta, double e) {
  double x = 0;
  double y = 0;
  double z = 0;
  double depth = 0.89*(7.7+ log(e) );
  double theta = 2*atan(exp(-1*eta));
  double r = 0;
  if( fabs(eta) > 1.479 )
    {
      double ecalZ = 315.4*fabs(eta)/eta;

      r = ecalZ / cos( 2*atan( exp( -1*eta ) ) ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  else
    {

      if(eta == 0) eta = 0.0001; // temporary measure to fix non when eta is exactly 0
      double rperp = 129.0;
      double zface = sqrt( cos( theta ) * cos( theta ) /
                            ( 1 - cos( theta ) * cos( theta ) ) *
                            rperp * rperp ) * fabs( eta ) / eta;
      r = sqrt( rperp * rperp + zface * zface ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  GlobalPoint pos(x,y,z);
  return pos;
}

double L1PiXTRKTreeProducer::etaToTheta(double eta)
{
  //  if(eta<0) return -2*atan(exp(eta));
  //  else return 2*atan(exp(-1*eta));
  return 2*atan(exp(-1*eta));
  //else return 2*atan(exp(-1*eta));

}

double L1PiXTRKTreeProducer::thetaToEta(double theta)
{
  //first bounds check theta to get into -pi/2 - pi/2 range
  while( fabs(theta) > L1PiXTRKTreeProducer::kPi/2.){
    if(theta>0) theta-=L1PiXTRKTreeProducer::kPi;
    else theta+=L1PiXTRKTreeProducer::kPi;
  }
  //now check sign
  if(theta<0) return log(tan(fabs(theta/2.)));
  else return -1.*log(tan(theta/2.));
}

double L1PiXTRKTreeProducer::detEtaFromEvnt(double evntEta,double z0){
  double thetaEvt = L1PiXTRKTreeProducer::etaToTheta(evntEta);
  double z = 129.4 / tan(thetaEvt); //129.4 is the average barrel radius
  double zTot = z+z0;

  if(fabs(zTot)<269){ //269 is an emperically derived number which means that < its likely in th barrel
    return zTot !=0 ? L1PiXTRKTreeProducer::thetaToEta(atan(129.4/zTot)) : 0.; //otherwise endcap time
  }
  double endcapZ = 319.2; //average z position of endcap
  if(evntEta<0) endcapZ*=-1;
  double rxy = tan(thetaEvt) * (endcapZ-z0);
  return L1PiXTRKTreeProducer::thetaToEta(atan(rxy/endcapZ));

}

// ------------ method called once each job just before starting event loop  ------------
void 
L1PiXTRKTreeProducer::beginJob(void)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1PiXTRKTreeProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void
L1PiXTRKTreeProducer::beginRun(edm::Run const& run, edm::EventSetup const& iSetup)
{  
   edm::ESHandle<HepPDT::ParticleDataTable> pdt;
   iSetup.getData(pdt);

   iSetup.get<IdealMagneticFieldRecord>().get(magField_);

}

unsigned int L1PiXTRKTreeProducer::getLayerNumber(const DetId& detid, const TrackerTopology* topo) {
    if (detid.det() == DetId::Tracker) {
        if (detid.subdetId() == PixelSubdetector::PixelBarrel) return (topo->pxbLayer(detid));
        else if (detid.subdetId() == PixelSubdetector::PixelEndcap) return (100 * topo->pxfSide(detid) + topo->pxfDisk(detid));
        else return 999;
    }
    return 999;
}



// From /RecoEgamma/EgammaMCTools/src/ElectronMCTruthFinder.cc
void L1PiXTRKTreeProducer::fill(const std::vector<SimTrack>& simTracks,
                                 const std::vector<SimVertex>& simVertices ) {
  //std::cout << "  ElectronMCTruthFinder::fill " << std::endl;

  unsigned nVtx = simVertices.size();
  unsigned nTks = simTracks.size();

  // Empty event, do nothin'
  if ( nVtx == 0 ) return;

  // create a map associating geant particle id and position in the 
  // event SimTrack vector
  for( unsigned it=0; it<nTks; ++it ) {
    geantToIndex_[ simTracks[it].trackId() ] = it;
    //std::cout << " ElectronMCTruthFinder::fill it " << it << " simTracks[it].trackId() " <<  simTracks[it].trackId() << std::endl;

  }
}

void
L1PiXTRKTreeProducer::doTrackMatching(const l1t::HGCalMulticluster& cluster, edm::Handle<L1TkTrackCollectionType> l1trackHandle)
{
  // track matching stuff
  // match to closes track up to delta R = 0.3
  // then match to the highest pt track < 0.3
  double min_track_dr = 999.;
  double max_track_pt = -999.;
  double max_track_pt_all_tracks = -999.;
  double max_track_pt_all_tracksEta = -999.;
  double max_track_pt_all_tracksPhi = -999.;
  double max_track_pt_all_tracksChi2 = -999.;
  double max_track_pt_all_chi2_cut = -999.;
  double max_track_pt_all_chi2_cutEta = -999.;
  double max_track_pt_all_chi2_cutPhi = -999.;
  double max_track_pt_all_chi2_cutChi2 = -999.;
  double matched_z = 999.;
  //edm::Ptr<TTTrack<Ref_PixelDigi_>> matched_track;
  bool foundMatchedTrack = false;
  edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> matched_track;
  if ( l1trackHandle.isValid() )
  {
     for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
     {
        //edm::Ptr<TTTrack<Ref_PixelDigi_>> ptr(l1trackHandle, track_index);
        edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> ptr(l1trackHandle, track_index);
        double pt = ptr->getMomentum().perp();

        // Don't consider tracks with pt < 2 for studies, might be increased to 3 later
        if (pt < 2.) continue;

        // Record the highest pt track per event
        // And, highest pt passing chi2 cut
        // to see if brem is an issues for mis-matched tracks
        double chi2 = ptr->getChi2();
        if (pt > max_track_pt_all_tracks) {
          max_track_pt_all_tracks = pt;
          max_track_pt_all_tracksEta = ptr->getMomentum().eta();
          max_track_pt_all_tracksPhi = ptr->getMomentum().phi();
          max_track_pt_all_tracksChi2 = chi2;}
        if (pt > max_track_pt_all_chi2_cut ) { // && chi2 < 100) {
          max_track_pt_all_chi2_cut = pt;
          max_track_pt_all_chi2_cutEta = ptr->getMomentum().eta();
          max_track_pt_all_chi2_cutPhi = ptr->getMomentum().phi();
          max_track_pt_all_chi2_cutChi2 = chi2;}

        // L1 Tracks are considered mis-measured if pt > 50
        // Therefore pt -> 50 if pt > 50
        // Only consider tracks if chi2 < 100
        //double dr = .1;

        GlobalPoint gp_cluster(cluster.position().x(), cluster.position().y(), cluster.position().z());

        double dr = L1TkElectronTrackMatchAlgo::deltaR(gp_cluster, ptr);
        //if (pt > 50.) pt = 50;
        // Choose closest track until dR < 0.3
        if ( dr < min_track_dr && min_track_dr > 0.3 ) // && chi2 < 100. )
        {
           min_track_dr = dr;
           max_track_pt = pt;
           matched_track = ptr;
           foundMatchedTrack = true;
           matched_z = ptr->getPOCA().z();
        }
        // If dR < 0.3, choose highest pt track
        else if ( dr < 0.3 && pt > max_track_pt ) // && chi2 < 100. )
        {
           min_track_dr = dr;
           max_track_pt = pt;
           matched_track = ptr;
           foundMatchedTrack = true;
           matched_z = ptr->getPOCA().z();
        }
     }
     float isoConeTrackCount = 0.; // matched track will be in deltaR cone
     float isoConePtSum = 0.;
     for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
     {
        //edm::Ptr<TTTrack<Ref_PixelDigi_>> ptr(l1trackHandle, track_index);
        edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> ptr(l1trackHandle, track_index);
        // don't double count the matched_track
        if ( ptr == matched_track ) {
          continue;
        }
        // Don't consider tracks with pt < 2 for studies, might be increased to 3 later
        double pt = ptr->getMomentum().perp();
        if (pt < 2.) continue;

        // Track Isolation cuts have been updated based on the recommendations here:
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1TrackTriggerObjects62X#Matching_L1EG_objects_with_L1Tra
        // Inner dR cone of 0.03, outer dR cone 0.2, momentum at least 2GeV
        // L1 Tracks are considered mis-measured if pt > 50
        // Therefore pt -> 50 if pt > 50
        // Only consider tracks if chi2 < 100
        // Only consider iso tracks from within dZ < 0.6
        double dr_2 = 99.;
        if (foundMatchedTrack) {
           dr_2 = reco::deltaR(ptr->getMomentum(), matched_track->getMomentum());}
        double chi2 = ptr->getChi2();
        //if (pt > 50.) pt = 50;
        double this_z = ptr->getPOCA().z();
        double deltaZ = abs(matched_z - this_z);
        if ( dr_2 < 0.2 && dr_2 > 0.03 && pt > 2. && deltaZ < 0.6 ) // chi2 < 100 && deltaZ < 0.6 )
        {
          isoConeTrackCount++;
          isoConePtSum += pt;
        }
     }

     hgcal_trackHighestPt.push_back(max_track_pt_all_tracks);
     hgcal_trackHighestPtEta.push_back(max_track_pt_all_tracksEta);
     hgcal_trackHighestPtPhi.push_back(max_track_pt_all_tracksPhi);
     hgcal_trackHighestPtChi2.push_back(max_track_pt_all_tracksChi2);

     hgcal_trackHighestPtCutChi2.push_back(max_track_pt_all_chi2_cut);
     hgcal_trackHighestPtCutChi2Eta.push_back(max_track_pt_all_chi2_cutEta);
     hgcal_trackHighestPtCutChi2Phi.push_back(max_track_pt_all_chi2_cutPhi);
     hgcal_trackHighestPtCutChi2Chi2.push_back(max_track_pt_all_chi2_cutChi2);

     hgcal_isTrackMatched.push_back(foundMatchedTrack);
     hgcal_isoConeNTrack.push_back(isoConeTrackCount);
     hgcal_isoConePtTrack.push_back(isoConePtSum);
     hgcal_trackmatchingdR.push_back(min_track_dr);

  }
}



void
L1PiXTRKTreeProducer::doTrackMatching(const l1slhc::L1EGCrystalCluster& cluster, edm::Handle<L1TkTrackCollectionType> l1trackHandle)
{
  // track matching stuff
  // match to closes track up to delta R = 0.3
  // then match to the highest pt track < 0.3
  double min_track_dr = 999.;
  double max_track_pt = -999.;
  double max_track_pt_all_tracks = -999.;
  double max_track_pt_all_tracksEta = -999.;
  double max_track_pt_all_tracksPhi = -999.;
  double max_track_pt_all_tracksChi2 = -999.;
  double max_track_pt_all_chi2_cut = -999.;
  double max_track_pt_all_chi2_cutEta = -999.;
  double max_track_pt_all_chi2_cutPhi = -999.;
  double max_track_pt_all_chi2_cutChi2 = -999.;
  double matched_z = 999.;
  //edm::Ptr<TTTrack<Ref_PixelDigi_>> matched_track;
  bool foundMatchedTrack = false;
  edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> matched_track;
  if ( l1trackHandle.isValid() )
  {
     for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
     {
        //edm::Ptr<TTTrack<Ref_PixelDigi_>> ptr(l1trackHandle, track_index);
        edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> ptr(l1trackHandle, track_index);
        double pt = ptr->getMomentum().perp();

        // Don't consider tracks with pt < 2 for studies, might be increased to 3 later
        if (pt < 2.) continue;

        // Record the highest pt track per event
        // And, highest pt passing chi2 cut
        // to see if brem is an issues for mis-matched tracks
        double chi2 = ptr->getChi2();
        if (pt > max_track_pt_all_tracks) {
          max_track_pt_all_tracks = pt;
          max_track_pt_all_tracksEta = ptr->getMomentum().eta();
          max_track_pt_all_tracksPhi = ptr->getMomentum().phi();
          max_track_pt_all_tracksChi2 = chi2;}
        if (pt > max_track_pt_all_chi2_cut ) { // && chi2 < 100) {
          max_track_pt_all_chi2_cut = pt;
          max_track_pt_all_chi2_cutEta = ptr->getMomentum().eta();
          max_track_pt_all_chi2_cutPhi = ptr->getMomentum().phi();
          max_track_pt_all_chi2_cutChi2 = chi2;}
        

        // L1 Tracks are considered mis-measured if pt > 50
        // Therefore pt -> 50 if pt > 50
        // Only consider tracks if chi2 < 100
        //double dr = .1;
        double dr = L1TkElectronTrackMatchAlgo::deltaR(L1TkElectronTrackMatchAlgo::calorimeterPosition(cluster.phi(), cluster.eta(), cluster.energy()), ptr);
        //if (pt > 50.) pt = 50;
        // Choose closest track until dR < 0.3
        if ( dr < min_track_dr && min_track_dr > 0.3 ) // && chi2 < 100. )
        {
           min_track_dr = dr;
           max_track_pt = pt;
           matched_track = ptr;
           foundMatchedTrack = true;
           matched_z = ptr->getPOCA().z();
        }
        // If dR < 0.3, choose highest pt track
        else if ( dr < 0.3 && pt > max_track_pt ) // && chi2 < 100. )
        {
           min_track_dr = dr;
           max_track_pt = pt;
           matched_track = ptr;
           foundMatchedTrack = true;
           matched_z = ptr->getPOCA().z();
        }
     }
     float isoConeTrackCount = 0.; // matched track will be in deltaR cone
     float isoConePtSum = 0.;
     for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
     {
        //edm::Ptr<TTTrack<Ref_PixelDigi_>> ptr(l1trackHandle, track_index);
        edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> ptr(l1trackHandle, track_index);
        // don't double count the matched_track
        if ( ptr == matched_track ) {
          continue;
        }
        // Don't consider tracks with pt < 2 for studies, might be increased to 3 later
        double pt = ptr->getMomentum().perp();
        if (pt < 2.) continue;

        // Track Isolation cuts have been updated based on the recommendations here:
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1TrackTriggerObjects62X#Matching_L1EG_objects_with_L1Tra
        // Inner dR cone of 0.03, outer dR cone 0.2, momentum at least 2GeV
        // L1 Tracks are considered mis-measured if pt > 50
        // Therefore pt -> 50 if pt > 50
        // Only consider tracks if chi2 < 100
        // Only consider iso tracks from within dZ < 0.6
        double dr_2 = 99.;
        if (foundMatchedTrack) {
           dr_2 = reco::deltaR(ptr->getMomentum(), matched_track->getMomentum());}
        double chi2 = ptr->getChi2();
        //if (pt > 50.) pt = 50;
        double this_z = ptr->getPOCA().z();
        double deltaZ = abs(matched_z - this_z);
        if ( dr_2 < 0.2 && dr_2 > 0.03 && pt > 2. && deltaZ < 0.6 ) // chi2 < 100 && deltaZ < 0.6 )
        {
          isoConeTrackCount++;
          isoConePtSum += pt;
        }
     }

     trackHighestPt.push_back(max_track_pt_all_tracks);
     trackHighestPtEta.push_back(max_track_pt_all_tracksEta);
     trackHighestPtPhi.push_back(max_track_pt_all_tracksPhi);
     trackHighestPtChi2.push_back(max_track_pt_all_tracksChi2);
     trackHighestPtCutChi2.push_back(max_track_pt_all_chi2_cut);
     trackHighestPtCutChi2Eta.push_back(max_track_pt_all_chi2_cutEta);
     trackHighestPtCutChi2Phi.push_back(max_track_pt_all_chi2_cutPhi);
     trackHighestPtCutChi2Chi2.push_back(max_track_pt_all_chi2_cutChi2);

     isTrackMatched.push_back(foundMatchedTrack);
     isoConeNTrack.push_back(isoConeTrackCount);
     isoConePtTrack.push_back(isoConePtSum);
     trackmatchingdR.push_back(min_track_dr);

  }
}

int L1PiXTRKTreeProducer::getLastCopy(int eleIndx, edm::Handle<reco::GenParticleCollection>& genParticles)
{
  int jda1 = -1;
  int jda2 = -1;
  const reco::GenParticle & genPart = (*genParticles)[eleIndx]; 
  int nrDa = genPart.numberOfDaughters();

  auto daughter = static_cast<const reco::GenParticle*>(genPart.daughter(0));
  // index of the first and the last daughter
  if(daughter){
    jda1 = std::distance(&(*genParticles)[0],daughter);
    jda2 = jda1 + nrDa-1;
  }

  // loop over daughters
  for(int i = jda1; i <= jda2; i++){
     const reco::GenParticle & part = (*genParticles)[i];
     if(part.pdgId() == genPart.pdgId()) {
       return getLastCopy(i, genParticles);
     } 
  }

 return eleIndx;
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1PiXTRKTreeProducer);
