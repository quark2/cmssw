
/** \Class GEMMaskReClusterizer
 *  \author J.C. Sanabria -- UniAndes, Bogota
 */

#include "GEMMaskReClusterizer.h"



GEMMaskReClusterizer::GEMMaskReClusterizer()
{

}


GEMMaskReClusterizer::~GEMMaskReClusterizer()
{

}


GEMClusterContainer GEMMaskReClusterizer::doAction(const GEMDetId& id,
                                                    GEMClusterContainer& initClusters,
                                                    const EtaPartitionMask& mask)
{

  GEMClusterContainer finClusters;
  GEMCluster prev;

  unsigned int j = 0;


  for (GEMClusterContainer::const_iterator i = initClusters.begin(); i != initClusters.end(); i++ ) {

    GEMCluster cl = *i;

    if ( i == initClusters.begin() ) {
      prev = cl;
      j++;
      if ( j == initClusters.size() ) {
	finClusters.insert(prev);
      }
      else if ( j < initClusters.size() ) {
	continue;
      }
    }


    if ( (prev.firstStrip()-cl.lastStrip()) == 2 and this->get(mask,cl.lastStrip()+1)
	 and cl.bx() == prev.bx() ) {

      GEMCluster merged(cl.firstStrip(),prev.lastStrip(),cl.bx());
      prev = merged;
      j++;
      if ( j == initClusters.size() ) {
	finClusters.insert(prev);
      }
    }

    else {
      j++;
      if ( j < initClusters.size() ) {
	finClusters.insert(prev);
	prev = cl;
      }
      if ( j == initClusters.size() ) {
	finClusters.insert(prev);
	finClusters.insert(cl);
      }
    }
  }

  return finClusters;

}



bool GEMMaskReClusterizer::get(const EtaPartitionMask& mask, int strip)
{
  return mask.test(strip-1);
}
