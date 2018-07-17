// -*- C++ -*-
//
// Package:    METExtractor/METAnalyzer
// Class:      METAnalyzer
// 
/**\class METAnalyzer METAnalyzer.cc METExtractor/METAnalyzer/plugins/METAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Byeonghak Ko
//         Created:  Mon, 16 Jul 2018 09:14:25 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/MET.h"

#include <TH1F.h>
#include <TFile.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class METAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit METAnalyzer(const edm::ParameterSet&);
    ~METAnalyzer();

    //static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

     // ----------member data ---------------------------
  private:
    edm::EDGetTokenT<std::vector<pat::MET>> tokenPatMET_;
    
    TH1F* m_h1Phi;
    TH1F* m_h1PhiXY;
    TH1F* m_h1Phi1XY;
    TH1F* m_h1Phi01XY;
    TH1F* m_h1Phi1SmearXY;
    TH1F* m_h1Phi01SmearXY;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
METAnalyzer::METAnalyzer(const edm::ParameterSet& iConfig)
  : tokenPatMET_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("src")))

{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  
  m_h1Phi = fs->make<TH1F>("phi_corr", "phi_corr", 40, -3.141592, 3.141593);
  m_h1PhiXY = fs->make<TH1F>("phi_corrXY", "phi_corrXY", 40, -3.141592, 3.141593);
  m_h1Phi1XY = fs->make<TH1F>("phi_corr1XY", "phi_corr1XY", 40, -3.141592, 3.141593);
  m_h1Phi01XY = fs->make<TH1F>("phi_corr01XY", "phi_corr01XY", 40, -3.141592, 3.141593);
  m_h1Phi1SmearXY = fs->make<TH1F>("phi_corr1SmearXY", "phi_corr1SmearXY", 40, -3.141592, 3.141593);
  m_h1Phi01SmearXY = fs->make<TH1F>("phi_corr01SmearXY", "phi_corr01SmearXY", 40, -3.141592, 3.141593);

}


METAnalyzer::~METAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
// 
// ------------ method called for each event  ------------
int g_nEvt = 0;
void
METAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  g_nEvt++;
   
  edm::Handle<std::vector<pat::MET>> srcPatMETCollection;
  iEvent.getByToken(tokenPatMET_, srcPatMETCollection);
  
  //std::cout << "MET (" << srcPatMETCollection->size() << ") : " 
  //  << ( *srcPatMETCollection )[ 0 ].phi() << std::endl;
  //if ( srcPatMETCollection->size() > 1 ) 
  //  printf("///////////////////////////!!!!///////////////// (%i)\n", g_nEvt);
  m_h1Phi->Fill(( *srcPatMETCollection )[ 0 ].phi());
  m_h1PhiXY->Fill(( *srcPatMETCollection )[ 0 ].corPhi(pat::MET::TypeXY));
  m_h1Phi1XY->Fill(( *srcPatMETCollection )[ 0 ].corPhi(pat::MET::Type1XY));
  m_h1Phi01XY->Fill(( *srcPatMETCollection )[ 0 ].corPhi(pat::MET::Type01XY));
  m_h1Phi1SmearXY->Fill(( *srcPatMETCollection )[ 0 ].corPhi(pat::MET::Type1SmearXY));
  m_h1Phi01SmearXY->Fill(( *srcPatMETCollection )[ 0 ].corPhi(pat::MET::Type01SmearXY));
}


// ------------ method called once each job just before starting event loop  ------------
void 
METAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*void
METAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}*/

//define this as a plug-in
DEFINE_FWK_MODULE(METAnalyzer);
