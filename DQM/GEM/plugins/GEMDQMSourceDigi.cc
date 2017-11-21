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

//#include "DataFormats/GEMRecHit/interface/GEMRecHit.h"
//#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMVfatStatusDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMGEBStatusDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMCStatusDigiCollection.h"


#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"


#include <string>

//----------------------------------------------------------------------------------------------------
 
class GEMDQMSourceDigi: public DQMEDAnalyzer
{
public:
  GEMDQMSourceDigi(const edm::ParameterSet& cfg);
  ~GEMDQMSourceDigi() override;
  
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

  edm::EDGetToken tagDigi;
  edm::EDGetToken tagError;
  edm::EDGetToken tagGEB;
  edm::EDGetToken tagAMC;

  const GEMGeometry* initGeometry(edm::EventSetup const & iSetup);
  int findVFAT(float min_, float max_, float x_, int roll_);
     
  const GEMGeometry* GEMGeometry_; 

  std::vector<GEMChamber> gemChambers;

  std::unordered_map<UInt_t,  MonitorElement*> Digi_Strip_vs_eta;
  //std::unordered_map<UInt_t,  MonitorElement*> h1B1010;
  //std::unordered_map<UInt_t,  MonitorElement*> h1B1100;
  //std::unordered_map<UInt_t,  MonitorElement*> h1B1110;
  //std::unordered_map<UInt_t,  MonitorElement*> h1Flag;
  
  MonitorElement *h1B1010All;
  MonitorElement *h1B1100All;
  MonitorElement *h1B1110All;
  
  MonitorElement *h1FlagAll;
  MonitorElement *h1CRCAll;
  
  MonitorElement *h1InputID;
  MonitorElement *h1Vwh;
  MonitorElement *h1Vwt;
  
  MonitorElement *h1GEBError;
  MonitorElement *h1GEBWarning;

  MonitorElement *GEMDAV;  
  MonitorElement *Tstate;  
  MonitorElement *GDcount; 
  MonitorElement *ChamT;   
  MonitorElement *OOSG;    


};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

int GEMDQMSourceDigi::findVFAT(float min_, float max_, float x_, int roll_) {
  float step = abs(max_-min_)/3.0;
  if ( x_ < (min(min_,max_)+step) ) { return 8 - roll_;}
  else if ( x_ < (min(min_,max_)+2.0*step) ) { return 16 - roll_;}
  else { return 24 - roll_;}
}

const GEMGeometry* GEMDQMSourceDigi::initGeometry(edm::EventSetup const & iSetup) {
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
GEMDQMSourceDigi::GEMDQMSourceDigi(const edm::ParameterSet& cfg)
{

  tagDigi = consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("digisInputLabel")); 
  tagError = consumes<GEMVfatStatusDigiCollection>(cfg.getParameter<edm::InputTag>("errorsInputLabel")); 
  tagGEB = consumes<GEMGEBStatusDigiCollection>(cfg.getParameter<edm::InputTag>("GEBInputLabel")); 
  tagAMC = consumes<GEMAMCStatusDigiCollection>(cfg.getParameter<edm::InputTag>("AMCInputLabel")); 

}

//----------------------------------------------------------------------------------------------------

GEMDQMSourceDigi::~GEMDQMSourceDigi()
{
}

//----------------------------------------------------------------------------------------------------

void GEMDQMSourceDigi::dqmBeginRun(edm::Run const &, edm::EventSetup const &)
{
}

//----------------------------------------------------------------------------------------------------

void GEMDQMSourceDigi::bookHistograms(DQMStore::IBooker &ibooker, edm::Run const &, edm::EventSetup const & iSetup)
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
  ibooker.setCurrentFolder("GEM/digi");
  for (auto ch : gemChambers){
    GEMDetId gid = ch.id();
    string hName_digi = "Digi_Strips_Gemini_"+to_string(gid.chamber())+"_la_"+to_string(gid.layer());
    string hTitle_digi = "Digi Strip Gemini ID : "+to_string(gid.chamber())+", layer : "+to_string(gid.layer());
    //     string hName_digi = "digi_"+to_string(gid.chamber());
    //     string hTitle_digi = "digi "+to_string(gid.chamber());
    Digi_Strip_vs_eta[ ch.id() ] = ibooker.book2D(hName_digi, hTitle_digi, 384, 0.5, 384.5, 8, 0.5,8.5);
    //string hNameErrors = "vfatErrors_"+to_string(gid.chamber())+"_la_"+to_string(gid.layer());
    //h1B1010[ ch.id() ] = ibooker.book1D(hNameErrors+"_b1010", hNameErrors+"_b1010", 15, 0x0 , 0xf);   
    //h1B1100[ ch.id() ] = ibooker.book1D(hNameErrors+"_b1100", hNameErrors+"_b1100", 15, 0x0 , 0xf);   
    //h1B1110[ ch.id() ] = ibooker.book1D(hNameErrors+"_b1110", hNameErrors+"_b1110", 15, 0x0 , 0xf);   
    //h1Flag[ ch.id() ] = ibooker.book1D(hNameErrors+"_flag", hNameErrors+"_flag", 15, 0x0 , 0xf);   
    //TH2F *hist_3 = Digi_Strip_vs_eta[ ch.id() ]->getTH2F();
    //hist_3->SetMarkerStyle(20);
    //hist_3->SetMarkerSize(0.5);
  }
  
  h1B1010All = ibooker.book1D("vfatErrors_all_b1010", "Control Bit 1010", 15, 0x0 , 0xf);   
  h1B1100All = ibooker.book1D("vfatErrors_all_b1100", "Control Bit 1100", 15, 0x0 , 0xf);   
  h1B1110All = ibooker.book1D("vfatErrors_all_b1110", "Control Bit 1110", 15, 0x0 , 0xf);   
  
  h1FlagAll = ibooker.book1D("vfatErrors_all_flag", "Control Flags", 15, 0x0 , 0xf);   
  h1CRCAll = ibooker.book1D("vfatErrors_all_CRC", "CRC Mismatches", 0xffff, -32768, 32768);   
  
  h1InputID = ibooker.book1D("GEB_InputID", "GEB GLIB input ID", 31,  0x0 , 0b11111);
  h1Vwh = ibooker.book1D("VFAT_Vwh", "VFAT word count", 4095,  0x0 , 0xfff);
  h1Vwt = ibooker.book1D("VFAT_Vwt", "VFAT word count", 4095,  0x0 , 0xfff);
  
  h1GEBError = ibooker.book1D("GEB_Errors", "GEB Critical Errors", 5, 0, 5);
  TH1F *histErr = h1GEBError->getTH1F();
  const char *error_flags[5] = {"Event Size Overflow", "L1AFIFO Full", "InFIFO Full", "Evt FIFO Full","InFIFO Underflow"};
  for (int i = 1; i<6; i++) histErr->GetXaxis()->SetBinLabel(i, error_flags[i-1]);
  h1GEBWarning = ibooker.book1D("GEB_Warnings", "GEB Warnings", 10,  0, 10);
  TH1F *histWar = h1GEBWarning->getTH1F();
  const char *warning_flags[10] = {"BX AMC-OH Mismatch", "BX AMC-VFAT Mismatch", "OOS AMC OH", "OOS AMC VFAT","No VFAT Marker","Event Size Warn", "L1AFIFO Near Full", "InFIFO Near Full", "EvtFIFO Near Full", "Stuck Data"};
  for (int i = 1; i<11; i++) histWar->GetXaxis()->SetBinLabel(i, warning_flags[i-1]);
  GEMDAV = ibooker.book1D("GEMDAV", "GEM DAV list", 24,  0, 24);
  Tstate     = ibooker.book1D("Tstate", "TTS state", 15,  0, 15);
  GDcount    = ibooker.book1D("GDcount", "GEM DAV count", 32,  0, 32);
  ChamT      = ibooker.book1D("ChamT", "Chamber Timeout", 24, 0, 24);
  OOSG       = ibooker.book1D("OOSG", "OOS GLIB", 1, 0, 1);
}

//----------------------------------------------------------------------------------------------------

void GEMDQMSourceDigi::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, 
                                            edm::EventSetup const& context) 
{
}

//----------------------------------------------------------------------------------------------------

void GEMDQMSourceDigi::analyze(edm::Event const& event, edm::EventSetup const& eventSetup)
{
  const GEMGeometry* GEMGeometry_  = initGeometry(eventSetup);
  if ( GEMGeometry_ == nullptr) return; 

  ////////////////
  ///// Digi /////
  ////////////////
  edm::Handle<GEMDigiCollection> gemDigis;
  edm::Handle<GEMVfatStatusDigiCollection> gemErrors;
  edm::Handle<GEMGEBStatusDigiCollection> gemGEB;
  edm::Handle<GEMAMCStatusDigiCollection> gemAMC;
  event.getByToken( this->tagDigi, gemDigis);
  event.getByToken( this->tagError, gemErrors);
  event.getByToken( this->tagGEB, gemGEB);
  event.getByToken( this->tagAMC, gemAMC);
  //   if (!gemDigis.isValid()){
  //   		edm::LogError("GEMDQMSourceDigi") << "GEM Digi is not valid.\n";
  //   		return;
  //   }
  for (auto ch : gemChambers){
    GEMDetId cId = ch.id();
     
    for(auto roll : ch.etaPartitions()){
      GEMDetId rId = roll->id();      
      const auto& digis_in_det = gemDigis->get(rId);
      for (auto d = digis_in_det.first; d != digis_in_det.second; ++d){
	Digi_Strip_vs_eta[ cId ]->Fill(d->strip(), rId.roll());
      }
      const auto& errors_in_det = gemErrors->get(rId);
      for(auto vfatError = errors_in_det.first; vfatError != errors_in_det.second; ++vfatError ){
        //h1B1010[ cId ]->Fill(vfatError->getB1010());
        //h1B1100[ cId ]->Fill(vfatError->getB1110());
        //h1B1110[ cId ]->Fill(vfatError->getB1110());
        //h1Flag[ cId ]->Fill(vfatError->getFlag());
        
        h1B1010All->Fill(vfatError->getB1010());
        h1B1100All->Fill(vfatError->getB1100());
        h1B1110All->Fill(vfatError->getB1110());
        h1FlagAll->Fill(vfatError->getFlag());
        h1CRCAll->Fill(vfatError->getCrc());
      }
    }
    const auto& GEB_in_det = gemGEB->get(cId);
    for(auto GEBStatus = GEB_in_det.first; GEBStatus != GEB_in_det.second; ++GEBStatus ){
      h1InputID->Fill(GEBStatus->getInputID());
      h1Vwh->Fill(GEBStatus->getVwh());
      h1Vwt->Fill(GEBStatus->getVwt());
      
      //h1GEBError->Fill(GEBStatus->getErrorC());
      for ( int bin = 0 ; bin < 9  ; bin++ ) 
        if ( ( ( GEBStatus->getErrorC() >> bin ) & 0x1 ) != 0 ) h1GEBWarning->Fill(bin);
      for ( int bin = 9 ; bin < 13 ; bin++ ) 
        if ( ( ( GEBStatus->getErrorC() >> bin ) & 0x1 ) != 0 ) h1GEBError->Fill(bin - 9);
      
      if ( ( GEBStatus->getInFu()   & 0x1 ) != 0 ) h1GEBError->Fill(9);
      if ( ( GEBStatus->getStuckd() & 0x1 ) != 0 ) h1GEBWarning->Fill(9);
    }
  }

  for (GEMAMCStatusDigiCollection::DigiRangeIterator amcIt = gemAMC->begin(); amcIt != gemAMC->end(); ++amcIt){
    const GEMAMCStatusDigiCollection::Range& range = (*amcIt).second;
    for ( auto amc = range.first; amc != range.second; ++amc ) {
    uint8_t binFired = 0;
    for (int bin = 0; bin < 24; bin++){
      binFired = ((amc->GEMDAV() >> bin) & 0x1);
      if (binFired) GEMDAV->Fill(bin);
      //binFired = ((amc->Bstatus() >> bin) & 0x1);
      //if (binFired) Bstatus->Fill(bin);
      binFired = ((amc->ChamT() >> bin) & 0x1);
      if (binFired) ChamT->Fill(bin);
    }
    Tstate->Fill(amc->Tstate());
    GDcount->Fill(amc->GDcount());
    OOSG->Fill(amc->OOSG());
  }
  }
}

//----------------------------------------------------------------------------------------------------

void GEMDQMSourceDigi::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& eSetup) 
{
}

//----------------------------------------------------------------------------------------------------

void GEMDQMSourceDigi::endRun(edm::Run const& run, edm::EventSetup const& eSetup)
{
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(GEMDQMSourceDigi);
