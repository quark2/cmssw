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

#include "DataFormats/GEMDigi/interface/GEMVfatStatusDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMGEBdataCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMCdataCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMC13EventCollection.h"

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
  void analyze(edm::Event const& e, edm::EventSetup const& eSetup) override;
  void endRun(edm::Run const& run, edm::EventSetup const& eSetup) override {};

private:
  int nVfat_ = 24;
  int cBit_ = 9;
  int eBit_ = 13;
  edm::EDGetToken tagVFAT_;
  edm::EDGetToken tagGEB_;
  edm::EDGetToken tagAMC_;
  edm::EDGetToken tagAMC13Event_;

  MonitorElement *h1_vfat_quality_;
  MonitorElement *h1_vfat_flag_;
  MonitorElement *h2_vfat_quality_;
  MonitorElement *h2_vfat_flag_;

  MonitorElement *h1_geb_inputStatus_;
  MonitorElement *h1_geb_vfatWordCnt_;
  MonitorElement *h1_geb_zeroSupWordsCnt_;
  MonitorElement *h1_geb_stuckData_;
  MonitorElement *h1_geb_inFIFOund_;
  
  MonitorElement *h1_amc_ttsState_;
  MonitorElement *h1_amc_davCnt_;
  MonitorElement *h1_amc_buffState_;
  MonitorElement *h1_amc_oosGlib_;
  MonitorElement *h1_amc_chTimeOut_;

  MonitorElement *h2_vfatPos_vs_channel_;
  MonitorElement *h1_amc13Event_bx_ok_;
  MonitorElement *h1_amcData_bx_ok_;

};

using namespace std;
using namespace edm;

GEMDQMStatusDigi::GEMDQMStatusDigi(const edm::ParameterSet& cfg)
{

  tagVFAT_ = consumes<GEMVfatStatusDigiCollection>(cfg.getParameter<edm::InputTag>("VFATInputLabel")); 
  tagGEB_ = consumes<GEMGEBdataCollection>(cfg.getParameter<edm::InputTag>("GEBInputLabel")); 
  tagAMC_ = consumes<GEMAMCdataCollection>(cfg.getParameter<edm::InputTag>("AMCInputLabel")); 
  tagAMC13Event_ = consumes<GEMAMC13EventCollection>(cfg.getParameter<edm::InputTag>("AMC13EventInputLabel"));

}

void GEMDQMStatusDigi::fillDescriptions(edm::ConfigurationDescriptions & descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("VFATInputLabel", edm::InputTag("muonGEMDigis", "vfatStatus")); 
  desc.add<edm::InputTag>("GEBInputLabel", edm::InputTag("muonGEMDigis", "gebStatus"));
  desc.add<edm::InputTag>("AMCInputLabel", edm::InputTag("muonGEMDigis", "AMCdata"));
  desc.add<edm::InputTag>("AMC13EventInputLabel", edm::InputTag("muonGEMDigis", "AMC13Event"));
  descriptions.add("GEMDQMStatusDigi", desc);  
}

void GEMDQMStatusDigi::bookHistograms(DQMStore::IBooker &ibooker, edm::Run const &, edm::EventSetup const & iSetup)
{

  ibooker.cd();
  ibooker.setCurrentFolder("GEM/StatusDigi");
  
  h1_vfat_quality_ = ibooker.book1D("vfat quality", "quality", 6, 0, 6);
  h1_vfat_flag_ = ibooker.book1D("vfat flag", "flag", 5, 0, 5);

  h2_vfat_quality_ = ibooker.book2D("vfat quality per geb", "quality", 6, 0 , 6, 36, 0, 36);
  h2_vfat_flag_ = ibooker.book2D("vfat flag per geb", "flag", 5, 0, 5, 36, 0, 36);

  h1_geb_inputStatus_ = ibooker.book1D("geb input status", "inputStatus", 10, 0, 10);
  h1_geb_vfatWordCnt_ = ibooker.book1D("geb no vfats", "nvfats", 25, 0, 25);
  h1_geb_zeroSupWordsCnt_ = ibooker.book1D("geb zeroSupWordsCnt", "zeroSupWordsCnt", 10, 0, 10);
  h1_geb_stuckData_ = ibooker.book1D("geb stuckData", "stuckData", 10, 0, 10);
  h1_geb_inFIFOund_ = ibooker.book1D("geb inFIFOund", "inFIFOund", 10, 0, 10);

  h1_amc_ttsState_ = ibooker.book1D("amc ttsState", "ttsState", 10, 0, 10);
  h1_amc_davCnt_ = ibooker.book1D("amc davCnt", "davCnt", 10, 0, 10);
  h1_amc_buffState_ = ibooker.book1D("amc buffState", "buffState", 10, 0, 10);
  h1_amc_oosGlib_ = ibooker.book1D("amc oosGlib", "oosGlib", 10, 0, 10);
  h1_amc_chTimeOut_ = ibooker.book1D("amc chTimeOut", "chTimeOut", 10, 0, 10);

  h2_vfatPos_vs_channel_ = ibooker.book2D("vfatPos vs Channel","vfatPos_vs_Channel;vfatPos;channel", 24,-0.5,23.5, 128,-0.5,127.5);
  h1_amc13Event_bx_ok_ = ibooker.book1D("amc13Event bx ok","amc13Event bx_header-bx_trailer;bx_header - bx_trailer;count",7,-3.5,3.5);
  h1_amcData_bx_ok_ = ibooker.book1D("amcData bx ok","amcData bx ok;amc13Ev_bx_header - amcData_bxId;count",7,-3.5,3.5);


}

void GEMDQMStatusDigi::analyze(edm::Event const& event, edm::EventSetup const& eventSetup)
{
  edm::Handle<GEMVfatStatusDigiCollection> gemVFAT;
  edm::Handle<GEMGEBdataCollection> gemGEB;
  edm::Handle<GEMAMCdataCollection> gemAMC;
  edm::Handle<GEMAMC13EventCollection> gemAMC13Event;
  event.getByToken( tagVFAT_, gemVFAT);
  event.getByToken( tagGEB_, gemGEB);
  event.getByToken( tagAMC_, gemAMC);
  event.getByToken( tagAMC13Event_, gemAMC13Event);

  for (GEMVfatStatusDigiCollection::DigiRangeIterator vfatIt = gemVFAT->begin(); vfatIt != gemVFAT->end(); ++vfatIt){
    GEMDetId gemid = (*vfatIt).first;
    float nIdx = gemid.chamber() + (gemid.layer()-1)/2.0;
    const GEMVfatStatusDigiCollection::Range& range = (*vfatIt).second;
    for ( auto vfatStat = range.first; vfatStat != range.second; ++vfatStat ) {
        
      h1_vfat_quality_->Fill(vfatStat->quality());
      h1_vfat_flag_->Fill(vfatStat->flag());
      h2_vfat_quality_->Fill(vfatStat->quality(), nIdx);
      h2_vfat_flag_->Fill(vfatStat->flag(), nIdx);
    }
  }

  for (GEMGEBdataCollection::DigiRangeIterator gebIt = gemGEB->begin(); gebIt != gemGEB->end(); ++gebIt){
    const GEMGEBdataCollection::Range& range = (*gebIt).second;    
    for ( auto GEBStatus = range.first; GEBStatus != range.second; ++GEBStatus ) {

      h1_geb_inputStatus_->Fill(GEBStatus->inputStatus());
      h1_geb_vfatWordCnt_->Fill(GEBStatus->vfatWordCnt()/3);
      h1_geb_zeroSupWordsCnt_->Fill(GEBStatus->zeroSupWordsCnt());
      h1_geb_stuckData_->Fill(GEBStatus->stuckData());
      h1_geb_inFIFOund_->Fill(GEBStatus->inFIFOund());

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

  for (GEMAMC13EventCollection::DigiRangeIterator amc13EvIt = gemAMC13Event->begin(); amc13EvIt != gemAMC13Event->end(); ++amc13EvIt){
    const GEMAMC13EventCollection::Range& range = (*amc13EvIt).second;
    //std::cout << "gemAMC13Event empty range? " << ((range.first==range.second)? "yes":"no") << "\n";
    for ( auto amc13 = range.first; amc13 != range.second; ++amc13 ) {

      long int bxH = (long int)(amc13->get_bxId());
      long int bxT = (long int)(amc13->get_bxIdT());

      h1_amc13Event_bx_ok_ ->Fill( bxH - bxT );

      const std::vector<gem::AMCdata>* amcData= amc13->getAMCpayloads();
      if (!amcData) std::cout << "amcData is null" << std::endl;
      //std::cout << "amcData.size=" << amcData->size() << "\n";
      for (unsigned int i=0; i<amcData->size(); i++) {
	long int bxId= (long int)(amcData->at(i).get_bxID());
	h1_amcData_bx_ok_ ->Fill( bxH - bxId );

	const std::vector<gem::GEBdata> *gebs = amcData->at(i).gebs();
	if (!gebs) std::cout << "gebs are null" << std::endl;
	//else std::cout << "gebs are ok size=" << gebs->size() << std::endl;
	for (auto gebIt= gebs->begin(); gebIt!=gebs->end(); gebIt++) {
	  const std::vector<gem::VFATdata> *vfats = gebIt->vFATs();
	  if (!vfats) std::cout << "vfats are null" << std::endl;
	  //else std::cout << "vfats are ok size=" << vfats->size() << "\n";
	  for (auto vfatIt= vfats->begin(); vfatIt!=vfats->end(); vfatIt++) {
	    //std::cout << "check vfat 0x" << std::hex << vfatIt->lsData() << " 0x" << vfatIt->msData() << std::dec << "\n";
	    for (uint64_t iChan=0; iChan<128; iChan++) {
	      if (vfatIt->chanIsON(iChan)) {
		//std::cout << "vfat @ " << vfatIt->get_pos() << " chanOn=" << iChan << "\n";
		h2_vfatPos_vs_channel_ ->Fill(float(vfatIt->get_pos()), float(iChan));
	      }
	    }
	  }
	}
      }
    }
  }

}

DEFINE_FWK_MODULE(GEMDQMStatusDigi);
