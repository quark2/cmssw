#ifndef CondFormats_GEMObjects_GEMDeadStrips_h
#define CondFormats_GEMObjects_GEMDeadStrips_h

#include "CondFormats/Serialization/interface/Serializable.h"
#include <vector>
#include <iostream>

class GEMDeadStrips {

 public:
  struct DeadItem {
    int rawId;
    int strip;  
    COND_SERIALIZABLE;
  };
  
  GEMDeadStrips(){}
  ~GEMDeadStrips(){}

  std::vector<DeadItem> const & getDeadVec() const {return deadVec_;}
  std::vector<DeadItem> & getDeadVec() {return deadVec_;}
  void setDeadVec(std::vector<DeadItem> &deadVec) {deadVec_ = deadVec;}
  void setDeadVec(std::vector<DeadItem> const &deadVec) {deadVec_ = deadVec;}

 private:
  std::vector<DeadItem> deadVec_;

  COND_SERIALIZABLE;
};
#endif
