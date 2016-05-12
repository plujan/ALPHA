#include "PileupAnalyzer.h"


PileupAnalyzer::PileupAnalyzer() {
    // PU reweighting
    LumiWeights=new edm::LumiReWeighting("data/MC_True.root", "data/Prod6.root", "S10", "pileup");
    
    std::cout << " - PileupAnalyzer initialized" << std::endl;
}

PileupAnalyzer::~PileupAnalyzer() {
    delete LumiWeights;
}


// ---------- PILEUP ----------

float PileupAnalyzer::GetPUWeight(const edm::Event& iEvent) {
  //  int nPT(0);
  //  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchool2012PileupReweighting
  //  if(!iEvent.isRealData()) {
  //    edm::Handle<std::vector<PileupSummaryInfo> > PUInfo;
  //    iEvent.getByLabel(edm::InputTag("PileupSummaryInfo"), PUInfo);
  //    for(std::vector<PileupSummaryInfo>::const_iterator pvi=PUInfo->begin(), pvn=PUInfo->end(); pvi!=pvn; ++pvi) {
  //      if(pvi->getBunchCrossing()==0) nPT=pvi->getTrueNumInteractions(); // getPU_NumInteractions();
  //    }
  //    return LumiWeights->weight( nPT );
  //  }
    return 1.;
}

