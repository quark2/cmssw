#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"


#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "EventFilter/GEMRawToDigi/interface/GEMVfatStatusDigiCollection.h"
#include "EventFilter/GEMRawToDigi/interface/GEMGEBdataCollection.h"
#include "EventFilter/GEMRawToDigi/interface/GEMAMCdataCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"

#include <string>

//----------------------------------------------------------------------------------------------------
 
class GEMDQMStatusDigi: public DQMEDAnalyzer
{
public:
  GEMDQMStatusDigi(const edm::ParameterSet& cfg);
  ~GEMDQMStatusDigi() override {};
  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions); 

protected:
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void bookHistogramsChamberPart(DQMStore::IBooker &, GEMDetId &);
  void bookHistogramsStationPart(DQMStore::IBooker &, GEMDetId &);
  
  void analyze(edm::Event const& e, edm::EventSetup const& eSetup) override;
  void endRun(edm::Run const& run, edm::EventSetup const& eSetup) override;

private:
  const GEMGeometry* initGeometry(edm::EventSetup const & iSetup);
  
  void AddLabel();
  
  void FillBits(MonitorElement *monitor, uint64_t unVal, int nNumBits);
  void FillBits(MonitorElement *monitor, uint64_t unVal, int nNumBits, int nY);
  
  const GEMGeometry* GEMGeometry_; 
  std::vector<GEMChamber> gemChambers_;
  
  int nStation_ = 2; // 3 for GE1/1, GE2/1, ME0
  int nLayer_ = 2;
  int nVfat_ = 24;
  
  int cBit_ = 9;
  int qVFATBit_ = 5;
  int fVFATBit_ = 4;
  int eBit_ = 13;
  int gebOtherBit_ = 2;
  
  edm::EDGetToken tagVFAT_;
  edm::EDGetToken tagGEB_;
  edm::EDGetToken tagAMC_;
  edm::EDGetToken tagDigi_;
  
  std::vector<GEMDetId> m_listLayers;

  MonitorElement *h1_vfat_qualityflag_;
  MonitorElement *h2_vfat_qualityflag_;
  
  std::unordered_map<UInt_t, MonitorElement*> listVFATQualityFlag_;
  std::unordered_map<UInt_t, MonitorElement*> listVFATBC_;
  std::unordered_map<UInt_t, MonitorElement*> listVFATEC_;

  std::unordered_map<UInt_t, MonitorElement*> listGEBInputStatus_;
  std::unordered_map<UInt_t, MonitorElement*> listGEBInputID_;
  std::unordered_map<UInt_t, MonitorElement*> listGEBVFATWordCnt_;
  std::unordered_map<UInt_t, MonitorElement*> listGEBVFATWordCntT_;
  std::unordered_map<UInt_t, MonitorElement*> listGEBZeroSupWordsCnt_;
  std::unordered_map<UInt_t, MonitorElement*> listGEBbcOH_;
  std::unordered_map<UInt_t, MonitorElement*> listGEBecOH_;
  std::unordered_map<UInt_t, MonitorElement*> listGEBOHCRC_;
  
  std::unordered_map<UInt_t, Bool_t> m_mapStatusFill;
  std::unordered_map<UInt_t, Bool_t> m_mapStatusErr;
  
  MonitorElement *h1_amc_ttsState_;
  MonitorElement *h1_amc_davCnt_;
  MonitorElement *h1_amc_buffState_;
  MonitorElement *h1_amc_oosGlib_;
  MonitorElement *h1_amc_chTimeOut_;
  
  MonitorElement *m_summaryReport;
  MonitorElement *h2SummaryStatus;

};

const GEMGeometry* GEMDQMStatusDigi::initGeometry(edm::EventSetup const & iSetup) {
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

using namespace std;
using namespace edm;

GEMDQMStatusDigi::GEMDQMStatusDigi(const edm::ParameterSet& cfg)
{

  tagVFAT_ = consumes<GEMVfatStatusDigiCollection>(cfg.getParameter<edm::InputTag>("VFATInputLabel")); 
  tagGEB_ = consumes<GEMGEBdataCollection>(cfg.getParameter<edm::InputTag>("GEBInputLabel")); 
  tagAMC_ = consumes<GEMAMCdataCollection>(cfg.getParameter<edm::InputTag>("AMCInputLabel")); 
  tagDigi_ = consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("digisInputLabel")); 
  
}

void GEMDQMStatusDigi::fillDescriptions(edm::ConfigurationDescriptions & descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("VFATInputLabel", edm::InputTag("muonGEMDigis", "vfatStatus")); 
  desc.add<edm::InputTag>("GEBInputLabel", edm::InputTag("muonGEMDigis", "gebStatus")); 
  desc.add<edm::InputTag>("AMCInputLabel", edm::InputTag("muonGEMDigis", "AMCdata")); 
  desc.add<edm::InputTag>("digisInputLabel", edm::InputTag("muonGEMDigis", "")); 
  descriptions.add("GEMDQMStatusDigi", desc);  
}


void GEMDQMStatusDigi::bookHistogramsChamberPart(DQMStore::IBooker &ibooker, GEMDetId &gid) {
  std::string hName, hTitle;
  
  UInt_t unBinPos;
  
  std::string strIdxName  = "Gemini_" + to_string(gid.chamber()) + "_GE" + 
    ( gid.region() > 0 ? "p" : "m" ) + to_string(gid.station()) + "_" + to_string(gid.layer());
  std::string strIdxTitle = "GEMINIm" + to_string(gid.chamber()) + " in GE" + 
    ( gid.region() > 0 ? "+" : "-" ) + to_string(gid.station()) + "/" + to_string(gid.layer());
  
  hName = "vfatStatus_QualityFlag_" + strIdxName;
  hTitle = "VFAT quality " + strIdxTitle;
  hTitle += ";VFAT;";
  listVFATQualityFlag_[ gid ] = ibooker.book2D(hName, hTitle, 24, 0, 24, 9, 0, 9);
  
  hName = "vfatStatus_BC_" + strIdxName;
  hTitle = "VFAT bunch crossing " + strIdxTitle;
  hTitle += ";Bunch crossing;VFAT";
  listVFATBC_[ gid ] = ibooker.book2D(hName, hTitle, 10, -5, 5, 24, 0, 24);
  
  hName = "vfatStatus_EC_" + strIdxName;
  hTitle = "VFAT event counter " + strIdxTitle;
  hTitle += ";Event counter;VFAT";
  listVFATEC_[ gid ] = ibooker.book2D(hName, hTitle, 256, 0, 256, 24, 0, 24);
  
  unBinPos = 1;
  listVFATQualityFlag_[ gid ]->setBinLabel(unBinPos++, "Good", 2);
  listVFATQualityFlag_[ gid ]->setBinLabel(unBinPos++, "CRC fail", 2);
  listVFATQualityFlag_[ gid ]->setBinLabel(unBinPos++, "b1010 fail", 2);
  listVFATQualityFlag_[ gid ]->setBinLabel(unBinPos++, "b1100 fail", 2);
  listVFATQualityFlag_[ gid ]->setBinLabel(unBinPos++, "b1110 fail", 2);
  listVFATQualityFlag_[ gid ]->setBinLabel(unBinPos++, "Hamming error", 2);
  listVFATQualityFlag_[ gid ]->setBinLabel(unBinPos++, "AFULL", 2);
  listVFATQualityFlag_[ gid ]->setBinLabel(unBinPos++, "SEUlogic", 2);
  listVFATQualityFlag_[ gid ]->setBinLabel(unBinPos++, "SUEI2C", 2);
  
  m_mapStatusFill[ gid ] = false;
  m_mapStatusErr[ gid ] = false;
}


void GEMDQMStatusDigi::bookHistogramsStationPart(DQMStore::IBooker &ibooker, GEMDetId &lid) {
  UInt_t unBinPos;
  
  Int_t  re = lid.region();
  UInt_t st = lid.station();
  UInt_t la = lid.layer();
  
  auto newbookGEB = [](DQMStore::IBooker &ibooker, std::string strName, std::string strTitle, std::string strAxis, 
                       int nLayer, int nStation, int re, int nBin, float fMin, float fMax)->MonitorElement *
  {
    strName  = strName  + "_st_" + (re>=0 ? "p" : "m") + std::to_string(nStation) + "_la_" + std::to_string(nLayer);
    strTitle = strTitle + ", station : " + (re>=0 ? "+" : "-") + std::to_string(nStation) + 
      ", layer : " + std::to_string(nLayer) + ";Chamber;" + strAxis;
    
    return ibooker.book2D(strName, strTitle, 36, 0, 36, nBin, fMin, fMax);
  };
  
  listGEBInputStatus_[ lid ]  = newbookGEB(ibooker, "geb_input_status", 
      "inputStatus",       "",                           la, st, re, 15, 0, 15);
  listGEBInputID_[ lid ]      = newbookGEB(ibooker, "geb_input_ID", 
      "inputID",           "Input ID",                   la, st, re, 32, 0, 32);
  listGEBVFATWordCnt_[ lid ]  = newbookGEB(ibooker, "geb_no_vfats", 
      "nvfats in header",  "Number of VFATs in header",  la, st, re, 25, 0, 25);
  listGEBVFATWordCntT_[ lid ] = newbookGEB(ibooker, "geb_no_vfatsT", 
      "nvfats in trailer", "Number of VFATs in trailer", la, st, re, 25, 0, 25);
  listGEBZeroSupWordsCnt_[ lid ] = newbookGEB(ibooker, "geb_zeroSupWordsCnt", 
      "zeroSupWordsCnt",   "Zero sup. words count",      la, st, re, 10, 0, 10);
  
  listGEBbcOH_[ lid ]  = newbookGEB(ibooker, "geb_bcOH",  "OH bunch crossing", 
      "OH bunch crossing", la, st, re, 3600, 0, 3600);
  listGEBecOH_[ lid ]  = newbookGEB(ibooker, "geb_ecOH",  "OH event coounter", 
      "OH event counter", la, st, re, 256, 0, 256);
  listGEBOHCRC_[ lid ] = newbookGEB(ibooker, "geb_OHCRC", "CRC of OH data", 
      "CRC of OH data", la, st, re, 65536, 0, 65536);
  
  unBinPos = 1;
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "BX mismatch GLIB OH", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "BX mismatch GLIB VFAT", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "OOS GLIB OH", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "OOS GLIB VFAT", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "No VFAT marker", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "Event size warn", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "L1AFIFO near full", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "InFIFO near full", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "EvtFIFO near full", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "Event size overflow", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "L1AFIFO full", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "InFIFO full", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "EvtFIFO full", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "Stuck data", 2);
  listGEBInputStatus_[ lid ]->setBinLabel(unBinPos++, "FIFO underflow", 2);
}


void GEMDQMStatusDigi::bookHistograms(DQMStore::IBooker &ibooker, edm::Run const &, edm::EventSetup const & iSetup)
{
  GEMGeometry_ = initGeometry(iSetup);
  if ( GEMGeometry_ == nullptr) return;
  
  const std::vector<const GEMSuperChamber*>& superChambers_ = GEMGeometry_->superChambers();   
  for ( auto sch : superChambers_ ) {
    int nLayer = sch->nChambers();
    for ( int l = 0 ; l < nLayer ; l++ ) gemChambers_.push_back(*sch->chamber(l + 1));
  }
  
  ibooker.cd();
  ibooker.setCurrentFolder("GEM/StatusDigi");
  
  h1_vfat_qualityflag_ = ibooker.book1D("vfat_quality_flag", "quality and flag", 9, 0, 9);
  h2_vfat_qualityflag_ = ibooker.book2D("vfat_quality_flag_per_geb", "quality and flag", 36, 0, 36, 9, 0, 9);
  
  m_listLayers.clear();
  
  for ( auto ch : gemChambers_ ) {
    GEMDetId gid = ch.id();
    //if ( gid.station() != 1 ) continue; // Now only GE11 is needed
    bookHistogramsChamberPart(ibooker, gid);
    
    GEMDetId layerID(gid.region(), gid.ring(), gid.station(), gid.layer(), 0, 0);
    Bool_t bOcc = false;
    
    for ( auto lid : m_listLayers ) {
      if ( lid == layerID ) {
        bOcc = true;
        break;
      }
    }
    
    if ( !bOcc ) m_listLayers.push_back(layerID);
  }
  
  auto lambdaLayer = [](GEMDetId a, GEMDetId b)->Bool_t {
    Int_t nA = a.region() * ( 20 * a.station() + a.layer() );
    Int_t nB = b.region() * ( 20 * b.station() + b.layer() );
    return nA > nB;
  };
  
  std::sort(m_listLayers.begin(), m_listLayers.end(), lambdaLayer);
  
  for ( auto lid : m_listLayers ) {
    bookHistogramsStationPart(ibooker, lid);
  }
  
  h1_amc_ttsState_  = ibooker.book1D("amc_ttsState",  "ttsState",  10, 0, 10);
  h1_amc_davCnt_    = ibooker.book1D("amc_davCnt",    "davCnt",    10, 0, 10);
  h1_amc_buffState_ = ibooker.book1D("amc_buffState", "buffState", 10, 0, 10);
  h1_amc_oosGlib_   = ibooker.book1D("amc_oosGlib",   "oosGlib",   10, 0, 10);
  h1_amc_chTimeOut_ = ibooker.book1D("amc_chTimeOut", "chTimeOut", 10, 0, 10);
  
  ibooker.cd();
  //ibooker.setCurrentFolder("GEM/Summary");
  ibooker.setCurrentFolder("GEM/EventInfo");
  
  // TODO: We need a study for setting the rule for this
  m_summaryReport = ibooker.bookFloat("reportSummary");
  m_summaryReport->Fill(1.0);
  
  //h2SummaryStatus = ibooker.book2D("summary_statusChamber", ";Chamber;", 36, 0, 36, 
  h2SummaryStatus = ibooker.book2D("reportSummaryMap", ";Chamber;", 36, 0, 36, 
      m_listLayers.size(), 0, m_listLayers.size());
  
  Int_t nIdxLayer = 0;
  for ( auto lid : m_listLayers ) {
    std::string strLabel = std::string("GE") + ( lid.region() > 0 ? "+" : "-" ) + std::to_string(lid.station())
      + "/" + std::to_string(lid.layer());
    h2SummaryStatus->setBinLabel(nIdxLayer + 1, strLabel, 2);
    
    nIdxLayer++;
  }
}


void GEMDQMStatusDigi::FillBits(MonitorElement *monitor, uint64_t unVal, int nNumBits) {
  int i = 0;
  uint64_t unFlag = 1;
  
  for ( ; i < nNumBits ; i++, unFlag <<= 1 ) {
    if ( ( unVal & unFlag ) != 0 ) {
      monitor->Fill(i);
    }
  }
}


void GEMDQMStatusDigi::FillBits(MonitorElement *monitor, uint64_t unVal, int nNumBits, int nX) {
  int i = 0;
  uint64_t unFlag = 1;
  
  for ( ; i < nNumBits ; i++, unFlag <<= 1 ) {
    if ( ( unVal & unFlag ) != 0 ) {
      monitor->Fill(nX, i);
    }
  }
}


void GEMDQMStatusDigi::analyze(edm::Event const& event, edm::EventSetup const& eventSetup)
{
  edm::Handle<GEMVfatStatusDigiCollection> gemVFAT;
  edm::Handle<GEMGEBdataCollection> gemGEB;
  edm::Handle<GEMAMCdataCollection> gemAMC;
  edm::Handle<GEMDigiCollection> gemDigis;
  
  event.getByToken( tagVFAT_, gemVFAT);
  event.getByToken( tagGEB_, gemGEB);
  event.getByToken( tagAMC_, gemAMC);
  event.getByToken( tagDigi_, gemDigis);

  for (GEMVfatStatusDigiCollection::DigiRangeIterator vfatIt = gemVFAT->begin(); vfatIt != gemVFAT->end(); ++vfatIt){
    GEMDetId gemid = (*vfatIt).first;
    GEMDetId gemchId = gemid.chamberId();
    
    float nIdx = gemid.chamber() + (gemid.layer()-1)/2.0;
    int nRoll = gemid.roll();
    const GEMVfatStatusDigiCollection::Range& range = (*vfatIt).second;
    
    for ( auto vfatStat = range.first; vfatStat != range.second; ++vfatStat ) {
      uint64_t unQFVFAT = vfatStat->quality() | ( vfatStat->flag() << qVFATBit_ );
      if ( ( unQFVFAT & ~0x1 ) == 0 )  {
        unQFVFAT |= 0x1; // If no error, then it should be 'Good'
      } else { // Error!!
        m_mapStatusErr[ gemchId ] = true;
      }
      
      FillBits(h1_vfat_qualityflag_, unQFVFAT, qVFATBit_ + fVFATBit_);
      FillBits(h2_vfat_qualityflag_, unQFVFAT, qVFATBit_ + fVFATBit_, nIdx);
      
      int nVFAT = ( 8 - nRoll ) + 8 * vfatStat->phi(); // vfatStat.phi()? or 2 - vfatStat.phi()?
      
      FillBits(listVFATQualityFlag_[ gemchId ], unQFVFAT, qVFATBit_ + fVFATBit_, nVFAT);
      //listVFATBC_[ gemchId ]->Fill(vfatStat->bc(), nVFAT);
      listVFATEC_[ gemchId ]->Fill(vfatStat->ec(), nVFAT);
    }
  }

  for (GEMGEBdataCollection::DigiRangeIterator gebIt = gemGEB->begin(); gebIt != gemGEB->end(); ++gebIt){
    GEMDetId gemid = (*gebIt).first;
    GEMDetId lid(gemid.region(), gemid.ring(), gemid.station(), gemid.layer(), 0, 0);
    UInt_t unCh = gemid.chamber() - 1;
    
    const GEMGEBdataCollection::Range& range = (*gebIt).second;    
    for ( auto GEBStatus = range.first; GEBStatus != range.second; ++GEBStatus ) {
      uint64_t unStatus = GEBStatus->inputStatus() | 
        ( GEBStatus->stuckData() << eBit_ ) | ( GEBStatus->inFIFOund() << (eBit_+1) );
      
      std::cout << gemid;
      printf("%06X\n", (unsigned int)unStatus);
      if ( unStatus != 0 ) { // Error!
        m_mapStatusErr[ gemid ] = true;
      }
      
      FillBits(listGEBInputStatus_[ lid ], unStatus, eBit_ + gebOtherBit_, unCh);
      
      listGEBInputID_[ lid ]->Fill(unCh, GEBStatus->inputID());
      listGEBVFATWordCnt_[ lid ]->Fill(unCh, GEBStatus->vfatWordCnt()/3);
      listGEBVFATWordCntT_[ lid ]->Fill(unCh, GEBStatus->vfatWordCntT()/3);
      listGEBZeroSupWordsCnt_[ lid ]->Fill(unCh, GEBStatus->zeroSupWordsCnt());
      
      listGEBbcOH_[ lid ]->Fill(unCh, GEBStatus->bcOH());
      listGEBecOH_[ lid ]->Fill(unCh, GEBStatus->ecOH());
      listGEBOHCRC_[ lid ]->Fill(unCh, GEBStatus->ohcrc());
    }
  }

  for (GEMAMCdataCollection::DigiRangeIterator amcIt = gemAMC->begin(); amcIt != gemAMC->end(); ++amcIt){
    const GEMAMCdataCollection::Range& range = (*amcIt).second;
    for ( auto amc = range.first; amc != range.second; ++amc ) {
      h1_amc_ttsState_->Fill(amc->ttsState());
      h1_amc_davCnt_->Fill(amc->davCnt());
      h1_amc_buffState_->Fill(amc->buffState());
      h1_amc_oosGlib_->Fill(amc->oosGlib());
      h1_amc_chTimeOut_->Fill(amc->chTimeOut());
    }
  }
  
  auto findVFAT = [](float min_, float max_, float x_, int roll_)->int {
    float step = max_/3;
    if ( x_ < (min_+step) ) { return 8 - roll_;}
    else if ( x_ < (min_+2*step) )  { return 16 - roll_;}
    else { return 24 - roll_;}
  };
  
  // Checking if there is a fire (data)
  for ( auto ch : gemChambers_ ) {
    GEMDetId cId = ch.id();
    Bool_t bIsHit = false;
    
    // Because every fired strip in a same VFAT shares a same bx, we keep bx from only one strip
    std::unordered_map<Int_t, Int_t> mapBXVFAT;
    
    for ( auto roll : ch.etaPartitions() ) {
      GEMDetId rId = roll->id();      
      const auto &digis_in_det = gemDigis->get(rId);
      
      for ( auto d = digis_in_det.first ; d != digis_in_det.second ; ++d ){
        bIsHit = true;
        mapBXVFAT[ findVFAT(1, roll->nstrips(), d->strip(), rId.roll()) ] = d->bx();
      }
      
      if ( bIsHit ) break;
    }
    
    for ( auto bx : mapBXVFAT ) listVFATBC_[ cId ]->Fill(bx.second, bx.first);
    
    if ( bIsHit ) { // Data occur!
      m_mapStatusFill[ cId ] = true;
    }
  }
}


void GEMDQMStatusDigi::endRun(edm::Run const& run, edm::EventSetup const& eSetup) {
  auto seekIdx = [](std::vector<GEMDetId> &listLayers, GEMDetId &id)->Int_t {
    for ( Int_t nIdx = 0 ; nIdx < (Int_t)listLayers.size() ; nIdx++ ) if ( id == listLayers[ nIdx ] ) return nIdx;
    return -1;
  };
  
  std::cout << "End run" << std::endl;
  std::unordered_map<UInt_t, Int_t> mapDone;
  
  for ( auto chamber : gemChambers_ ) {
    auto gid = chamber.id();
    
    UInt_t unVal = 0; // No data, no error
    if      ( m_mapStatusErr[ gid ]  ) unVal = 1; // Error! (no matter there are data or not)
    else if ( m_mapStatusFill[ gid ] ) unVal = 2; // Data occur, with no error
    std::cout << gid << ": " << unVal << std::endl;
    
    GEMDetId layerId(gid.region(), gid.ring(), gid.station(), gid.layer(), 0, 0);
    Int_t nIdxLayer = seekIdx(m_listLayers, layerId);
    
    h2SummaryStatus->setBinContent(gid.chamber(), nIdxLayer + 1, unVal);
    
    // The below is for under/overflow of BX histograms
    Int_t nNbinX = listVFATBC_[ gid ]->getNbinsX();
    Int_t nNbinY = listVFATBC_[ gid ]->getNbinsY();
    
    for ( Int_t i = 0 ; i < nNbinY ; i++ ) {
      listVFATBC_[ gid ]->setBinContent(1, i, 
          listVFATBC_[ gid ]->getBinContent(0, i) + listVFATBC_[ gid ]->getBinContent(1, i));
      listVFATBC_[ gid ]->setBinContent(nNbinX, i, 
          listVFATBC_[ gid ]->getBinContent(nNbinX, i) + listVFATBC_[ gid ]->getBinContent(nNbinX + 1 , i));
    }
  }
}

DEFINE_FWK_MODULE(GEMDQMStatusDigi);
