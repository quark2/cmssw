#include "CondTools/GEM/interface/GEMMaskedStripSourceHandler.h"
#include "CondCore/CondDB/interface/ConnectionPool.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "RelationalAccess/ISessionProxy.h"
#include "RelationalAccess/ITransaction.h"
#include "RelationalAccess/ISchema.h"
#include "RelationalAccess/ITable.h"
#include "RelationalAccess/IQuery.h"
#include "RelationalAccess/ICursor.h"
#include "CoralBase/AttributeList.h"
#include "CoralBase/Attribute.h"
#include "CoralBase/AttributeSpecification.h"
#include <TString.h>

#include <fstream>
#include <cstdlib>
#include <vector>

#include <DataFormats/MuonDetId/interface/GEMDetId.h>

popcon::GEMMaskedStripSourceHandler::GEMMaskedStripSourceHandler( const edm::ParameterSet& ps ):
  m_name( ps.getUntrackedParameter<std::string>( "name", "GEMMaskedStripSourceHandler" ) ),
  m_dummy( ps.getUntrackedParameter<int>( "WriteDummy", 0 ) ),
  m_validate( ps.getUntrackedParameter<int>( "Validate", 1 ) ),
  m_connect( ps.getParameter<std::string>( "connect" ) ),
  m_connectionPset( ps.getParameter<edm::ParameterSet>( "DBParameters" ) )
{
}

popcon::GEMMaskedStripSourceHandler::~GEMMaskedStripSourceHandler()
{
}

void popcon::GEMMaskedStripSourceHandler::getNewObjects()
{
  
  edm::LogInfo( "GEMMaskedStripSourceHandler" ) << "[" << "GEMMaskedStripSourceHandler::" << __func__ << "]:" << m_name << ": "
                                         << "BEGIN" << std::endl;
  
  edm::Service<cond::service::PoolDBOutputService> mydbservice;
  
  // first check what is already there in offline DB
  Ref payload;
  if(m_validate==1 && tagInfo().size>0) {
    payload = lastPayload();
  }
  
  // now construct new cabling map from online DB
  // FIXME: use boost::ptime
  /*time_t rawtime;
  time(&rawtime); //time since January 1, 1970
  tm * ptm = gmtime(&rawtime);//GMT time
  char buffer[20];
  strftime(buffer,20,"%d/%m/%Y_%H:%M:%S",ptm);
  std::string eMap_version( buffer );*/
  stripsMap =  new GEMMaskedStrips();
  
  std::string baseCMS = std::string(getenv("CMSSW_BASE"))+std::string("/src/CondTools/GEM/data/");  
  std::vector<std::string> mapfiles;
  mapfiles.push_back("GEMMaskVec.dat");
  
  // VFAT Postion Map 
  std::string filename(baseCMS+mapfiles[0]);
  std::ifstream maptype(filename.c_str());
  
  std::string field, line;
  while(std::getline(maptype, line)){
    GEMMaskedStrips::MaskItem itemNew;
    
    std::stringstream ssline(line);
    getline( ssline, field, ' ' );
    std::stringstream IdRoll_Dead(field);
    getline( ssline, field, ' ' );
    std::stringstream IdStrip_Dead(field);
    
    IdRoll_Dead >> itemNew.rawId;
    IdStrip_Dead >> itemNew.strip;
    
    //std::cout<<" Sector="<<sec<<" z_direction="<<z_dir<<" ieta="<<ieta<<" iphi="<<iphi<<" depth="<<dep<<" vfat position="<<vfat_pos<<" vfat address = " << vfat_add << " vfat type = " << vfat_type << " amc ID = " << amc_ID << " geb ID = " << geb_ID << std::endl;
    std::cout << " ROLL ID=" << itemNew.rawId << " STRIP=" << itemNew.strip << std::endl;
    
    stripsMap->getMaskVec().push_back(itemNew);
  }
    
  cond::Time_t snc = mydbservice->currentTime();  
  // look for recent changes
  int difference=1;
  if (difference==1) {
    m_to_transfer.push_back(std::make_pair((GEMMaskedStrips*)stripsMap,snc));
  }
}

// // additional work (I added these two functions: ConnectOnlineDB and DisconnectOnlineDB)
void popcon::GEMMaskedStripSourceHandler::ConnectOnlineDB( const std::string& connect, const edm::ParameterSet& connectionPset )
{
  cond::persistency::ConnectionPool connection;
  connection.setParameters( connectionPset );
  connection.configure();
  session = connection.createSession( connect,true );
}

void popcon::GEMMaskedStripSourceHandler::DisconnectOnlineDB()
{
  session.close();
}
