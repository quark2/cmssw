/*
 *  \author Julia Yarba
 */

#include <ostream>

#include "IOMC/ParticleGuns/interface/FlatRandomPtGunProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "CLHEP/Random/RandFlat.h"

using namespace edm;
using namespace std;

FlatRandomPtGunProducer::FlatRandomPtGunProducer(const ParameterSet& pset) : 
   BaseFlatGunProducer(pset)
{


   ParameterSet defpset ;
   ParameterSet pgun_params = 
      pset.getParameter<ParameterSet>("PGunParameters") ;
  
   fMinPt = pgun_params.getParameter<double>("MinPt");
   fMaxPt = pgun_params.getParameter<double>("MaxPt");
  
  produces<HepMCProduct>("unsmeared");
  produces<GenEventInfoProduct>();
}

FlatRandomPtGunProducer::~FlatRandomPtGunProducer()
{
   // no need to cleanup GenEvent memory - done in HepMCProduct
}


int myIsMuonPassScint(double dVx, double dVy, double dVz, double dPx, double dPy, double dPz) {
  // To test the drop-down of efficiency at edges, we set the cut looser
  //double dL = -1200.0;
  //double dR =  1200.0;
  double dL = -300.0;
  double dR =  300.0;
  double dB =  -850.0;
  double dT =   950.0;
  //double dT =  -500.0;
  
  double dYLower = -114.85;
  double dYUpper = 1540.15;
  
  double dTLower = ( dYLower - dVy ) / dPy;
  double dTUpper = ( dYUpper - dVy ) / dPy;
  
  double dXLower = dVx + dTLower * dPx;
  double dZLower = dVz + dTLower * dPz;
  
  if ( !( dL <= dXLower && dXLower <= dR && dB <= dZLower && dZLower <= dT ) ) {
    return 0;
  }
  
  double dXUpper = dVx + dTUpper * dPx;
  double dZUpper = dVz + dTUpper * dPz;
  
  if ( !( dL <= dXUpper && dXUpper <= dR && dB <= dZUpper && dZUpper <= dT ) ) {
    return 0;
  }
  
  return 1;
}

void FlatRandomPtGunProducer::produce(Event &e, const EventSetup& es) 
{
   edm::Service<edm::RandomNumberGenerator> rng;
   CLHEP::HepRandomEngine* engine = &rng->getEngine(e.streamID());

   if ( fVerbosity > 0 )
   {
      cout << " FlatRandomPtGunProducer : Begin New Event Generation" << endl ; 
   }
   // event loop (well, another step in it...)
          
   // no need to clean up GenEvent memory - done in HepMCProduct
   // 
   
   // here re-create fEvt (memory)
   //
   fEvt = new HepMC::GenEvent() ;
   
   // now actualy, cook up the event from PDGTable and gun parameters
   //
   // 1st, primary vertex
   //
   //HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(0.,0.,0.));
   /*double dVx = CLHEP::RandFlat::shoot(engine, -1200.0, 1200.0) ;
   double dVy = -100.0;
   double dVz = CLHEP::RandFlat::shoot(engine, -650.0, 750.0) ;
   HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(dVx,dVy,dVz));*/
   double dVx;
   double dVy = 1550.0;
   double dVz;
   HepMC::GenVertex* Vtx = NULL;

   // loop over particles
   //
   int barcode = 1 ;
   for (unsigned int ip=0; ip<fPartIDs.size(); ++ip)
   {
       double px, py, pz, mom;
       int j = 0;
       
       // To avoid meeting a muon which does not pass the detectors apparently
       while ( 1 ) {
         //dVx = CLHEP::RandFlat::shoot(engine, -2300.0, 2300.0) ;
         dVx = CLHEP::RandFlat::shoot(engine, -300.0, 300.0) ;
         dVz = CLHEP::RandFlat::shoot(engine, -1000.0, 1000.0) ;
         //dVz = CLHEP::RandFlat::shoot(engine, -850.0, -450.0) ;
         
         double pt     = CLHEP::RandFlat::shoot(engine, fMinPt, fMaxPt) ;
         //double eta    = CLHEP::RandFlat::shoot(engine, fMinEta, fMaxEta) ;
         double phi    = CLHEP::RandFlat::shoot(engine, fMinPhi, fMaxPhi) ;
         //double theta  = 2.*atan(exp(-eta)) ;
         double theta  = CLHEP::RandFlat::shoot(engine, 0.0, 3.141592 / 2.0);
         //double theta  = acos(sqrt(CLHEP::RandFlat::shoot(engine, 0.0, 1.0))) ;
         //double theta  = acos(CLHEP::RandFlat::shoot(engine, 0.0, 1.0)) ;
         //double theta  = 0.0;
         //mom    = pt/sin(theta) ;
         mom    = pt;
         /*double px     = pt*cos(phi) ;
         double py     = pt*sin(phi) ;
         double pz     = mom*cos(theta) ;*/
         px     = mom*sin(theta)*cos(phi) ;
         py     = -mom*cos(theta) ;
         pz     = mom*sin(theta)*sin(phi) ;
         //px = 0.0;
         //py = mom;
         //pz = 0.0;
         
         if ( 1 == 1 || myIsMuonPassScint(dVx, dVy, dVz, px, py, pz) != 0 ) break;
         
         if ( j >= 10000 ) break;
         j++;
       }
       int PartID = fPartIDs[ip] ;
       const HepPDT::ParticleData* 
          PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID))) ;
       double mass   = PData->mass().value() ;
       Vtx = new HepMC::GenVertex(HepMC::FourVector(dVx,dVy,dVz));
       //px = 0;
       //py = pt;
       //pz = 0;
       double energy2= mom*mom + mass*mass ;
       double energy = sqrt(energy2) ; 
       HepMC::FourVector p(px,py,pz,energy) ;
       HepMC::GenParticle* Part = 
           new HepMC::GenParticle(p,PartID,1);
       Part->suggest_barcode( barcode ) ;
       barcode++ ;
       Vtx->add_particle_out(Part);

       if ( fAddAntiParticle )
       {
          HepMC::FourVector ap(-px,-py,-pz,energy) ;
	  int APartID = -PartID ;
	  if ( PartID == 22 || PartID == 23 )
	  {
	     APartID = PartID ;
	  }	  
	  HepMC::GenParticle* APart =
	     new HepMC::GenParticle(ap,APartID,1);
	  APart->suggest_barcode( barcode ) ;
	  barcode++ ;
	  Vtx->add_particle_out(APart) ;
       }

   }

   fEvt->add_vertex(Vtx) ;
   fEvt->set_event_number(e.id().event()) ;
   fEvt->set_signal_process_id(20) ; 
        
   if ( fVerbosity > 0 )
   {
      fEvt->print() ;  
   }

   unique_ptr<HepMCProduct> BProduct(new HepMCProduct()) ;
   BProduct->addHepMCData( fEvt );
   e.put(std::move(BProduct), "unsmeared");

   unique_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
   e.put(std::move(genEventInfo));
    
   if ( fVerbosity > 0 )
   {
      // for testing purpose only
      // fEvt->print() ; // prints empty info after it's made into edm::Event
      cout << " FlatRandomPtGunProducer : Event Generation Done " << endl;
   }
}
//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(FlatRandomPtGunProducer);
