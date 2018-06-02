
/** \Class GEMMaskReClusterizer
 *  \author J.C. Sanabria -- UniAndes, Bogota
 */

#include "GEMMaskReClusterizer.h"


GEMClusterContainer GEMMaskReClusterizer::doAction(const GEMDetId& id,
                                                   GEMClusterContainer& initClusters,
                                                   const EtaPartitionMask& mask) const
{
  GEMClusterContainer finClusters;
  if ( initClusters.empty() ) return finClusters;

  GEMCluster prev = *initClusters.begin();
  for ( auto cl = std::next(initClusters.begin()); cl != initClusters.end(); ++cl ) {
    // Merge this cluster if it is adjacent by 1 masked strip
    // Note that the GEMClusterContainer collection is sorted in DECREASING ORDER of strip #
    // So the prev. cluster is placed after the current cluster (check the < operator of GEMCluster carefully)
    if ( (prev.firstStrip()-cl->lastStrip()) == 2 and
         this->get(mask, cl->lastStrip()+1) and prev.bx() == cl->bx() ) {
      GEMCluster merged(cl->firstStrip(), prev.lastStrip(), cl->bx());
      prev = merged;
    }
    else {
      finClusters.insert(prev);
      prev = *cl;
    }
  }

  // Finalize by putting the last cluster to the collection
  finClusters.insert(prev);

  return finClusters;

}

bool GEMMaskReClusterizer::get(const EtaPartitionMask& mask, int strip) const
{
  return mask.test(strip-1);
}

/*GEMMaskReClusterizer::GEMMaskReClusterizer()
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


    if ( ((prev.firstStrip()-cl.lastStrip()) == 2 && this->get(mask,(cl.lastStrip()+1)))
	 && (cl.bx() == prev.bx()) ) {

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



int GEMMaskReClusterizer::get(const EtaPartitionMask& mask, int strip)
{

  if ( mask.test(strip-1) ) return 1;
  else return 0;

}*/
