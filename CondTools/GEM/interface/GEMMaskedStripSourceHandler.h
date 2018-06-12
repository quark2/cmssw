#ifndef GEMEMAPSOURCEHANDLER
#define GEMEMAPSOURCEHANDLER

#include "CondCore/PopCon/interface/PopConSourceHandler.h"
#include "CondFormats/GEMObjects/interface/GEMMaskedStrips.h"
#include "CondFormats/GEMObjects/interface/GEMDeadStrips.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>

namespace popcon
{
  class GEMMaskedStripSourceHandler : public popcon::PopConSourceHandler<GEMMaskedStrips>
  {
      
  public:
      
    GEMMaskedStripSourceHandler( const edm::ParameterSet& ps );
    ~GEMMaskedStripSourceHandler();
    void getNewObjects();
    void ConnectOnlineDB( const std::string& connect, const edm::ParameterSet& connectionPset ); // additional work
    void DisconnectOnlineDB(); // additional work
    void readGEMEMap();
    std::string id() const { return m_name; }
      
  private:
      
    GEMMaskedStrips * stripsMap;
    GEMDeadStrips * stripsMap2;
    cond::persistency::Session session; // additional work
    std::string m_name;
    int m_dummy;
    int m_validate;
    std::string m_connect;
    std::string m_authpath;
    edm::ParameterSet m_connectionPset;
  };
}
#endif
