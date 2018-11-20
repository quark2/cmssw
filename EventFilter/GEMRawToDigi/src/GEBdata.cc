#include "EventFilter/GEMRawToDigi/interface/GEBdata.h"
#include <string>
#include <sstream>

std::string gem::GEBdata::getChamberHeader_str() const {
  std::stringstream ss;
  ss << "ChamberHeader: unused1=" << (uint32_t)(ch_.unused1) << " inputStatus=" << (uint32_t)(ch_.inputStatus) << " [flag_noVFATmarker=" << flag_noVFATmarker() << "]" << " vfatWordCnt=" << (uint32_t)(ch_.vfatWordCnt) << " inputID=" << (uint32_t)(ch_.inputID) << " zeroSupWordsCnt=" << (uint32_t)(ch_.zeroSupWordsCnt) << " unused2=" << (uint32_t)(ch_.unused2);
  return ss.str();
}

std::string gem::GEBdata::getChamberTrailer_str() const {
  std::stringstream ss;
  ss << "ChamberTrailer: unused=" << (uint32_t)(ct_.unused) << " inUfw=" << (uint32_t)(ct_.inUfw) << " stuckData=" << (uint32_t)(ct_.stuckData) << " inFIFOund=" << (uint32_t)(ct_.inFIFOund) << " vfatWordCntT=" << (uint32_t)(ct_.vfatWordCntT) << " ohcrc=" << (uint32_t)(ct_.ohcrc);
  return ss.str();
}
