#include "GEMClusterizer.h"
#include "GEMCluster.h"
#include "GEMClusterContainer.h"


GEMClusterContainer GEMClusterizer::doAction(const GEMDigiCollection::Range& digiRange)
{
  GEMClusterContainer initialCluster, finalCluster;
  // Return empty container for null input
  if ( std::distance(digiRange.second, digiRange.first) == 0 ) return finalCluster;

  // Start from single digi recHits
  for ( auto digi = digiRange.first; digi != digiRange.second; ++digi ) {
    GEMCluster cl(digi->strip(), digi->strip(), digi->bx());
    //if ( digi->hasTime() ) cl.addTime(digi->time());
    //if ( digi->hasY() ) cl.addY(digi->coordinateY());
    initialCluster.insert(cl);
  }
  if ( initialCluster.empty() ) return finalCluster; // Confirm the collection is valid

  // Start from the first initial cluster
  GEMCluster prev = *initialCluster.begin();

  // Loop over the remaining digis
  // Note that the last one remains as open in this loop
  for ( auto cl = std::next(initialCluster.begin()); cl != initialCluster.end(); ++cl ) {
    if ( prev.isAdjacent(*cl) ) {
      // Merged digi to the previous one
      prev.merge(*cl);
    }
    else {
      // Close the previous cluster and start new cluster
      finalCluster.insert(prev);
      prev = *cl;
    }
  }

  // Finalize by adding the last cluster
  finalCluster.insert(prev);

  return finalCluster;
}

 
/*GEMClusterContainer
GEMClusterizer::doAction(const GEMDigiCollection::Range& digiRange){
  GEMClusterContainer cls;
  for (GEMDigiCollection::const_iterator digi = digiRange.first;
       digi != digiRange.second;
       digi++) {
    GEMCluster cl(digi->strip(),digi->strip(),digi->bx());
    cls.insert(cl);
  }
  GEMClusterContainer clsNew =this->doActualAction(cls);
  return clsNew;
}

GEMClusterContainer
GEMClusterizer::doActualAction(GEMClusterContainer& initialclusters){
  
  GEMClusterContainer finalCluster;
  GEMCluster prev;

  unsigned int j = 0;
  for(GEMClusterContainer::const_iterator i=initialclusters.begin();
      i != initialclusters.end(); i++){
    GEMCluster cl = *i;

    if(i==initialclusters.begin()){
      prev = cl;
      j++;
      if(j == initialclusters.size()){
	finalCluster.insert(prev);
      }
      else if(j < initialclusters.size()){
	continue;
      }
    }

    if(prev.isAdjacent(cl)) {
      prev.merge(cl);
      j++;
      if(j == initialclusters.size()){
	finalCluster.insert(prev);
      }
    }
    else {
      j++;
      if(j < initialclusters.size()){
	finalCluster.insert(prev);
	prev = cl;
      }
      if(j == initialclusters.size()){
	finalCluster.insert(prev);
	finalCluster.insert(cl);
      }
    }
  }

  return finalCluster;
}*/
 

