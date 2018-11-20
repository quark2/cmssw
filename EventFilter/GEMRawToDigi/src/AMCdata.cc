#include "EventFilter/GEMRawToDigi/interface/AMCdata.h"
#include <string>
#include <sstream>

std::string gem::AMCdata::getAMCheader1_str() const {
  std::stringstream ss;
  ss << "AMCheader1: reserved=" << (uint32_t)(amch1_.reserved) << " AMCnum=" << (uint32_t)(amch1_.AMCnum) << " l1AID=" << (uint32_t)(amch1_.l1AID) << " bxID=" << (uint32_t)(amch1_.bxID) << " dataLen=" << (uint32_t)(amch1_.dataLength);
  return ss.str();
}

std::string gem::AMCdata::getAMCheader2_str() const {
  std::stringstream ss;
  ss << "AMCheader2: boardID=" << (uint32_t)(amch2_.boardID) << " orbitNum=" << (uint32_t)(amch2_.orbitNum) << " param3=" << (uint32_t)(amch2_.param3) << " param2=" << (uint32_t)(amch2_.param2) << " param1=" << (uint32_t)(amch2_.param1) << " runType=" << (uint32_t)(amch2_.runType) << " formatVer=" << (uint32_t)(amch2_.formatVer);
  return ss.str();
}

std::string gem::AMCdata::getAMCtrailer_str() const {
  std::stringstream ss;
  ss << "AMCtrailer: dataLenT=" << (uint32_t)(amct_.dataLengthT) << " cb0=" << (uint32_t)(amct_.cb0) << " l1AIDT=" << (uint32_t)(amct_.l1AIDT) << " crc=" << (uint32_t)(amct_.crc);
  return ss.str();
}

std::string gem::AMCdata::getGEMeventHeader_str() const {
  std::stringstream ss;
  ss << "GEMeventHeader: ttsState=" << (uint32_t)(eh_.ttsState) << " reserved=" << (uint32_t)(eh_.reserved) << " davCnt=" << (uint32_t)(eh_.davCnt) << " buffState=" << (uint32_t)(eh_.buffState) << " davList=" << (uint32_t)(eh_.davList);
  return ss.str();
}

std::string gem::AMCdata::getGEMeventTrailer_str() const {
  std::stringstream ss;
  ss << "GEMeventTrailer: unused1=" << (uint32_t)(et_.unused1) << " ctrl5=" << (uint32_t)(et_.ctrl5) << " unused2=" << (uint32_t)(et_.unused2) << " oosGlib=" << (uint32_t)(et_.oosGlib) << " chTimeOut=" << (uint32_t)(et_.chTimeOut);
  return ss.str();
}
