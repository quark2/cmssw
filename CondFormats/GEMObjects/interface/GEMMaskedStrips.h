#ifndef CondFormats_GEMObjects_GEMMaskedStrips_h
#define CondFormats_GEMObjects_GEMMaskedStrips_h

#include "CondFormats/Serialization/interface/Serializable.h"
#include <vector>
#include <iostream>

class GEMMaskedStrips {

 public:
  struct MaskItem {
    int rawId;
    int strip;
    COND_SERIALIZABLE;
  };
  
  GEMMaskedStrips(){}
  ~GEMMaskedStrips(){}

  std::vector<MaskItem> const & getMaskVec() const {return maskVec_;}
  std::vector<MaskItem> & getMaskVec() {return maskVec_;}
  void setMaskVec(std::vector<MaskItem> &maskVec) {maskVec_ = maskVec;}
  void setMaskVec(std::vector<MaskItem> const &maskVec) {maskVec_ = maskVec;}

 private:
  std::vector<MaskItem> maskVec_;

  COND_SERIALIZABLE;
};
#endif
