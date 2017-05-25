/****************************************************************************

****************************************************************************/

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/GEMRecHit/interface/GEMRecHit.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"


#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"


#include <string>

//----------------------------------------------------------------------------------------------------
 
class GEMDQMSource: public DQMEDAnalyzer
{
  public:
    GEMDQMSource(const edm::ParameterSet& cfg);
    virtual ~GEMDQMSource();
  
  protected:
    void dqmBeginRun(edm::Run const &, edm::EventSetup const &) override;
    void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
    void analyze(edm::Event const& e, edm::EventSetup const& eSetup) override;
    void beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup) override;
    void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup) override;
    void endRun(edm::Run const& run, edm::EventSetup const& eSetup) override;

  private:
    unsigned int verbosity;
   
    int nCh;


    //edm::EDGetTokenT<  > gemDigi;
    edm::EDGetToken tagRecHit;
    edm::EDGetToken tagDigi;

    const GEMGeometry* initGeometry(edm::EventSetup const & iSetup);
    int findVFAT(float min_, float max_, float x_, int roll_);
     
    const GEMGeometry* GEMGeometry_; 

    std::vector<GEMChamber> gemChambers;

    MonitorElement* testPhi;
    MonitorElement* testEta;
    std::unordered_map<UInt_t,  MonitorElement*> recHitME;
    std::unordered_map<UInt_t,  MonitorElement*> digiME;

     

};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

int GEMDQMSource::findVFAT(float min_, float max_, float x_, int roll_) {
  float step = abs(max_-min_)/3.0;
  if ( x_ < (min(min_,max_)+step) ) { return 8 - roll_;}
  else if ( x_ < (min(min_,max_)+2.0*step) ) { return 16 - roll_;}
  else { return 24 - roll_;}
}

const GEMGeometry* GEMDQMSource::initGeometry(edm::EventSetup const & iSetup) {
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


//----------------------------------------------------------------------------------------------------
GEMDQMSource::GEMDQMSource(const edm::ParameterSet& cfg)
{

  tagRecHit = consumes<GEMRecHitCollection>(cfg.getParameter<edm::InputTag>("recHitsInputLabel")); 
  //tagDigi = consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("digisInputLabel")); 

}

//----------------------------------------------------------------------------------------------------

GEMDQMSource::~GEMDQMSource()
{
}

//----------------------------------------------------------------------------------------------------

void GEMDQMSource::dqmBeginRun(edm::Run const &, edm::EventSetup const &)
{
}

//----------------------------------------------------------------------------------------------------

void GEMDQMSource::bookHistograms(DQMStore::IBooker &ibooker, edm::Run const &, edm::EventSetup const & iSetup)
{

  GEMGeometry_ = initGeometry(iSetup);
  if ( GEMGeometry_ == nullptr) return ;  

  const std::vector<const GEMSuperChamber*>& superChambers_ = GEMGeometry_->superChambers();   
  for (auto sch : superChambers_){
    int n_lay = sch->nChambers();
    for (int l=0;l<n_lay;l++){
      gemChambers.push_back(*sch->chamber(l+1));
    }
  }
  nCh = gemChambers.size();


  ibooker.cd();
  ibooker.setCurrentFolder("GEM");
  testPhi = ibooker.book1D("testPhi", "testPhi", 360, -3.141592, 3.141592);
  testEta = ibooker.book1D("testEta", "testEta", 10, 1.5, 2.5);
  ibooker.setCurrentFolder("GEM/recHit");
  for (auto ch : gemChambers){
    GEMDetId gid = ch.id();
    string hName = "recHit_"+to_string(gid.chamber());
    string hTitle = "recHit "+to_string(gid.chamber());
    recHitME[ ch.id() ] = ibooker.book1D(hName, hTitle, 24,0,24);
  }
  /*
  ibooker.setCurrentFolder("GEM/digi");
  for (auto ch : gemChambers){
    GEMDetId gid = ch.id();
    string hName = "digi_"+to_string(gid.chamber());
    string hTitle = "digi "+to_string(gid.chamber());
    digiME[ ch.id() ] = ibooker.book1D(hName, hTitle, 24,0,24);
  }
  */


}

//----------------------------------------------------------------------------------------------------

void GEMDQMSource::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, 
                                            edm::EventSetup const& context) 
{
}

//----------------------------------------------------------------------------------------------------

void GEMDQMSource::analyze(edm::Event const& event, edm::EventSetup const& eventSetup)
{
  const GEMGeometry* GEMGeometry_  = initGeometry(eventSetup);
  if ( GEMGeometry_ == nullptr) return; 

  edm::Handle<GEMRecHitCollection> gemRecHits;
  event.getByToken( this->tagRecHit, gemRecHits);
  if (!gemRecHits.isValid()) {
    edm::LogError("GEMDQMSource") << "GEM recHit is not valid.\n";
    return ;
  }  
  for (GEMRecHitCollection::const_iterator recHit = gemRecHits->begin(); recHit != gemRecHits->end(); ++recHit){
    LocalPoint recHitLP = recHit->localPosition();
    GlobalPoint recHitGP = GEMGeometry_->idToDet((*recHit).gemId())->surface().toGlobal(recHitLP);
    Float_t rhPhi = recHitGP.phi();
    Float_t rhEta = recHitGP.eta();
    testPhi->Fill(rhPhi);
    testEta->Fill(rhEta);
  }
  for (auto ch : gemChambers){
    GEMDetId cId = ch.id();
    for(auto roll : ch.etaPartitions()){
      GEMDetId rId = roll->id();
      const auto& recHitsRange = gemRecHits->get(rId); 
      auto gemRecHit = recHitsRange.first;
      for ( auto hit = gemRecHit; hit != recHitsRange.second; ++hit ) {
        int nVfat = findVFAT(1.0, 385.0, hit->firstClusterStrip()+0.5*hit->clusterSize(), rId.roll());
        recHitME[ cId ]->Fill(nVfat);
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------

void GEMDQMSource::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& eSetup) 
{
}

//----------------------------------------------------------------------------------------------------

void GEMDQMSource::endRun(edm::Run const& run, edm::EventSetup const& eSetup)
{
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(GEMDQMSource);
