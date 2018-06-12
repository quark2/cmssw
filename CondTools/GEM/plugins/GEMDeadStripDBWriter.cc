#include "CondCore/PopCon/interface/PopConAnalyzer.h"
#include "CondTools/GEM/interface/GEMDeadStripSourceHandler.h"
#include "FWCore/Framework/interface/MakerMacros.h"

typedef popcon::PopConAnalyzer<popcon::GEMDeadStripSourceHandler> GEMDeadStripDBWriter;
//define this as a plug-in
DEFINE_FWK_MODULE(GEMDeadStripDBWriter);
