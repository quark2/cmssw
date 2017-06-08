#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//DQM services
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <memory>
#include <vector>

#include "TString.h"


using namespace std;

class GEMDQMHarvester: public edm::EDAnalyzer
{
  
 public:

  explicit GEMDQMHarvester(const edm::ParameterSet&);
  virtual ~GEMDQMHarvester(){};
  virtual void beginJob(){return;};
  virtual void endJob(){return;};  
  virtual void analyze(const edm::Event&, const edm::EventSetup&){return;};
  virtual void beginRun(const edm::Run&, const edm::EventSetup&){return;};
  virtual void endRun(const edm::Run&, const edm::EventSetup&){return;};

  
private:
  std::string fName;
  int verbosity;
  DQMStore *dbe;

};


GEMDQMHarvester::GEMDQMHarvester(const edm::ParameterSet& ps)
{
  fName = ps.getUntrackedParameter<std::string>("Name");

  //dbe_path_ = std::string("GEMDQM/");
  //outputFile_ = ps.getUntrackedParameter<std::string>("outputFile", "myfile.root");
}


//void GEMDQMHarvestor::dqmEndJob(DQMStore::IBooker & ibooker, DQMStore::IGetter &ig )
//{
  //ig.setCurrentFolder(dbe_path_.c_str());

//}
DEFINE_FWK_MODULE(GEMDQMHarvester);
