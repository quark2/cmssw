#ifndef EventFilter_GEMRawToDigi_AMC13Event_h
#define EventFilter_GEMRawToDigi_AMC13Event_h
#include <vector>
#include <string>
#include "AMCdata.h"

namespace gem {
  
  union CDFHeader {
    uint64_t word;
    struct {
      uint64_t fov       : 8;  // not used
      uint64_t sourceId  : 12; // FED number assigned by CDAQ
      uint64_t bxId      : 12; // Bunch crossing 0...3563
      uint64_t lv1Id     : 24; // Level 1 ID (hardware event number)
      uint64_t eventType : 4;  // Event Type (1 for normal, 2 for calibration)
      uint64_t cb5       : 4;  // 0x5
    };
  };
  union AMC13Header {
    uint64_t word;
    struct {
      uint64_t cb0       : 4;  // 0x0
      uint64_t orbitN    : 32; // Orbit Number
      uint64_t reserved0 : 16; // reserved
      uint64_t nAMC      : 4;  // Number of AMCs following (0 to 12)
      uint64_t calType   : 4;  // Calibration event type
      uint64_t uFov      : 4;  // Format version: 0x1
    };
  };

  union AMC13DataDef {
    AMC13DataDef(uint64_t d) : word(d) {}

    uint64_t word;
    struct {
      uint64_t boardId    : 16; // board id
      uint64_t amcNr      : 4;  // amcNr
      uint64_t blkSeqNo   : 8;  // block sequence no
      uint64_t cb0        : 4;  // 0x0
      uint64_t dataSize   : 24; // amc payload size
      uint64_t endBits8   : 8;  // c,v,p,e,s,m,l,?
    };
  };

  union AMC13Trailer {
    uint64_t word;
    struct {
      uint64_t bxIdT     : 12; // bx id
      uint64_t lv1IdT    : 8;  // level 1 id
      uint64_t blkN      : 8;  // block number
      uint64_t cb0       : 4;  // 0x0
      uint64_t crc32     : 32; // Overall CRC (first 32 bits)
    };
  };
  union CDFTrailer {
    uint64_t word;
    struct {
      uint64_t trdd      : 4;  // $,$,t,r
      uint64_t tts       : 4;  // tts (first 4 bits)
      uint64_t evtStat   : 4;  // event status
      uint64_t cfxx      : 4;  // x,x,f,c
      uint64_t crcCDF    : 16; // CDF crc (first 16 bits)
      uint64_t evtLength : 24; // event length
      uint64_t evtType   : 4;  // event type
      uint64_t cbA       : 4;  // 0xA (first 4 bits)
    };
  };

  class AMC13Event
  {
    
  public:
    AMC13Event() {}
    ~AMC13Event() {amcHeaders_.clear(); amcs_.clear();}

    void setCDFHeader(uint64_t word) { cdfh_.word = word;}
    void setCDFHeader(uint8_t Evt_ty, uint32_t LV1_id, uint16_t BX_id, uint16_t Source_id);
    uint64_t getCDFHeader() const { return cdfh_.word;}

    void setAMC13Header(uint64_t word) { amc13h_.word = word;}
    void setAMC13Header(uint8_t CalTyp, uint8_t nAMC, uint32_t OrN);
    uint64_t getAMC13Header() const { return amc13h_.word;}

    void setAMC13Trailer(uint64_t word) { amc13t_.word = word;}
    void setAMC13Trailer(uint8_t Blk_NoT, uint8_t LV1_idT, uint16_t BX_idT);
    uint64_t getAMC13Trailer() const { return amc13t_.word;}

    void setCDFTrailer(uint64_t word) { cdft_.word = word;}
    void setCDFTrailer(uint32_t EvtLength);
    uint64_t getCDFTrailer() const { return cdft_.word;}

    uint16_t bxId() const {return cdfh_.bxId;}
    uint32_t lv1Id() const {return cdfh_.lv1Id;}
    uint16_t sourceId() const {return cdfh_.sourceId;}

    uint32_t get_fov() const {return (uint32_t)cdfh_.fov;}
    uint32_t get_sourceId() const {return (uint32_t)cdfh_.sourceId;}
    uint32_t get_bxId() const {return (uint32_t)cdfh_.bxId;}
    uint32_t get_lv1Id() const { return (uint32_t)cdfh_.lv1Id;}
    uint32_t get_eventType() const {return (uint32_t)cdfh_.eventType;}
    uint32_t get_cb5() const {return (uint32_t)cdfh_.cb5;}

    uint8_t  nAMC() const {return amc13h_.nAMC;}

    uint32_t get_cb0() const {return (uint32_t)amc13h_.cb0;}
    uint32_t get_orbitN() const {return (uint32_t)amc13h_.orbitN;}
    uint32_t get_reserved0() const {return (uint32_t)amc13h_.reserved0;}
    uint32_t get_nAMC() const {return (uint32_t)amc13h_.nAMC;}
    uint32_t get_calType() const {return (uint32_t)amc13h_.calType;}
    uint32_t get_uFov() const {return (uint32_t)amc13h_.uFov;}

    const std::vector<uint64_t> * getAMCheaders() const {return &amcHeaders_;}
    void addAMCheader(uint64_t word);
    void addAMCheader(uint32_t AMC_size, uint8_t Blk_No, uint8_t AMC_No, uint16_t BoardID);

    uint32_t get_boardId(int i) const { AMC13DataDef d(amcHeaders_[i]); return (uint32_t)d.boardId; }
    uint32_t get_amcNr(int i) const { AMC13DataDef d(amcHeaders_[i]); return (uint32_t)d.amcNr; }
    uint32_t get_blkSeqNo(int i) const { AMC13DataDef d(amcHeaders_[i]); return (uint32_t)d.blkSeqNo; }
    uint32_t get_cb0(int i) const { AMC13DataDef d(amcHeaders_[i]); return (uint32_t)d.cb0; }
    uint32_t get_dataSize(int i) const { AMC13DataDef d(amcHeaders_[i]); return (uint32_t)d.dataSize; }
    uint32_t get_endBits8(int i) const { AMC13DataDef d(amcHeaders_[i]); return (uint32_t)d.endBits8; }

    uint32_t get_bxIdT() const {return (uint32_t)(amc13t_.bxIdT);}
    uint32_t get_lv1IdT() const {return (uint32_t)(amc13t_.lv1IdT);}
    uint32_t get_blkSeqNoT() const {return (uint32_t)(amc13t_.blkN);}
    uint32_t get_cb0T() const {return(uint32_t)(amc13t_.cb0);}

    uint32_t get_trddT() const {return (uint32_t)(cdft_.trdd);}
    uint32_t get_evtTypeT() const {return (uint32_t)(cdft_.evtType);}
    uint32_t get_evtLenT() const {return (uint32_t)(cdft_.evtLength);}
    uint32_t get_cbAT() const {return (uint32_t)(cdft_.cbA);}

    const std::vector<AMCdata> * getAMCpayloads() const {return &amcs_;}   
    void addAMCpayload(const AMCdata& a) {amcs_.push_back(a);}

    std::string getCDFHeader_str() const;
    std::string getCDFTrailer_str() const;
    std::string getAMC13Header_str() const;
    std::string getAMC13Trailer_str() const;
    std::string getAMC13DataDef_str(int i) const;
    
  private:
    // AMC headers
    std::vector<uint64_t> amcHeaders_;
    // AMCs payload
    std::vector<AMCdata> amcs_;

    CDFHeader cdfh_;
    AMC13Header amc13h_;
    AMC13Trailer amc13t_;
    CDFTrailer cdft_;
      
  };
}
#endif
