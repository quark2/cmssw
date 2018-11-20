#include <cstdint>
#include "EventFilter/GEMRawToDigi/interface/AMC13Event.h"
#include <sstream>

using namespace gem;

void AMC13Event::setCDFHeader(uint8_t Evt_ty, uint32_t LV1_id, uint16_t BX_id, uint16_t Source_id)
{
  cdfh_.cb5 = 0x5;
  cdfh_.eventType = Evt_ty;
  cdfh_.lv1Id = LV1_id;
  cdfh_.bxId = BX_id;
  cdfh_.sourceId = Source_id;
  cdfh_.fov = 0; // not used
}

void AMC13Event::setAMC13Header(uint8_t CalTyp, uint8_t nAMC, uint32_t OrN)
{
  amc13h_.uFov = 0; // format version # FIX ME
  amc13h_.calType = CalTyp;
  amc13h_.nAMC = nAMC;
  amc13h_.reserved0 = 0;
  amc13h_.orbitN = OrN;
  amc13h_.cb0 = 0x0;
}

void AMC13Event::setAMC13Trailer(uint8_t Blk_NoT, uint8_t LV1_idT, uint16_t BX_idT)
{
  amc13t_.crc32 = 0; // FIX ME
  amc13t_.blkN = Blk_NoT;
  amc13t_.lv1IdT = LV1_idT;
  amc13t_.bxIdT = BX_idT;
}

void AMC13Event::setCDFTrailer(uint32_t EvtLength)
{
  cdft_.cbA = 0xA;
  cdft_.evtType = 0; // FIX ME
  cdft_.evtLength = EvtLength;
  cdft_.crcCDF = 0; // FIX ME
  cdft_.cfxx = 0; // FIX ME
  cdft_.evtStat = 0; // event status # FIX ME
  cdft_.tts = 0;
  cdft_.trdd = 0;
}

void AMC13Event::addAMCheader(uint64_t word)
{
  amcHeaders_.push_back(word);
}

// FIX ME: not verified wrt to AMC13DataDef
void AMC13Event::addAMCheader(uint32_t AMC_size, uint8_t Blk_No, uint8_t AMC_No, uint16_t BoardID)
{
  // AMC Header word 
  // 55 - 32  | 27 - 20 | 19 - 16 | 15 - 0  |
  // AMC_size | Blk_No  | AMC_No  | BoardID |  
  uint64_t word =
    (static_cast<uint64_t>(AMC_size & 0x00ffffff) << 32) |
    (static_cast<uint64_t>(Blk_No & 0xff) << 20) |
    (static_cast<uint64_t>(AMC_No & 0x0f) << 16) |
    (static_cast<uint64_t>(BoardID & 0xffff));
  amcHeaders_.push_back(word);
}

std::string AMC13Event::getCDFHeader_str() const
{
  std::stringstream ss;
  ss << "CDFHeader: fov=" << (uint32_t)(cdfh_.fov) << " sourceId=" << (uint32_t)(cdfh_.sourceId) << " bxId=" << (uint32_t)(cdfh_.bxId) << " lv1Id=" << (uint32_t)(cdfh_.lv1Id) << " eventType=" << (uint32_t)(cdfh_.eventType) << " cb5=" << (uint32_t)(cdfh_.cb5);
  return ss.str();
}

std::string AMC13Event::getCDFTrailer_str() const
{
  std::stringstream ss;
  ss << "CDFTrailer: trdd=" << (uint32_t)(cdft_.trdd) << " tts=" << (uint32_t)(cdft_.tts) << " evtStat=" << (uint32_t)(cdft_.evtStat) << " cfxx=" << (uint32_t)(cdft_.cfxx) << " crcCDF=" << (uint32_t)(cdft_.crcCDF) << " evtLength=" << (uint32_t)(cdft_.evtLength) << " evtType=" << (uint32_t)(cdft_.evtType) << " cbA=" << (uint32_t)(cdft_.cbA);
  return ss.str();
}

std::string AMC13Event::getAMC13Header_str() const
{
  std::stringstream ss;
  ss << "AMC13Header: cb0=" << (uint32_t)(amc13h_.cb0) << " orbitN=" << (uint32_t)(amc13h_.orbitN) << " reserved0=" << (uint32_t)(amc13h_.reserved0) << " nAMC=" << (uint32_t)(amc13h_.nAMC) << " calType=" << (uint32_t)(amc13h_.calType) << " uFov=" << (uint32_t)(amc13h_.uFov);
  return ss.str();
}

std::string AMC13Event::getAMC13Trailer_str() const
{
  std::stringstream ss;
  ss << "AMC13Trailer: bxIdT=" << (uint32_t)(amc13t_.bxIdT) << " lv1IdT=" << (uint32_t)(amc13t_.lv1IdT) << " blkN=" << (uint32_t)(amc13t_.blkN) << " crc32=" << (uint32_t)(amc13t_.crc32);
  return ss.str();
}

std::string AMC13Event::getAMC13DataDef_str(int i) const
{
  AMC13DataDef d = amcHeaders_.at(i);
  std::stringstream ss;
  ss << "AMC13DataDef: boardId=" << (uint32_t)(d.boardId) << " amcNr=" << (uint32_t)(d.amcNr) << " blkSeqNo=" << (uint32_t)(d.blkSeqNo) << " cb0=" << (uint32_t)(d.cb0) << " dataSize=" << (uint32_t)(d.dataSize) << " endBits8" << (uint32_t)(d.endBits8);
  return ss.str();
}
