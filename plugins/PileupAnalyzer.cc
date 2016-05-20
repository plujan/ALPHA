#include "PileupAnalyzer.h"


PileupAnalyzer::PileupAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    PUToken(CColl.consumes<std::vector<PileupSummaryInfo> >(PSet.getParameter<edm::InputTag>("pileup"))),
    DataFileName(PSet.getParameter<std::string>("dataFileName")),
    MCFileName(PSet.getParameter<std::string>("mcFileName")),
    DataName(PSet.getParameter<std::string>("dataName")),
    MCName(PSet.getParameter<std::string>("mcName"))
{   
    std::cout << " --- PileupAnalyzer initialization ---" << std::endl;
    std::cout << "  pileup MC file    :\t" << MCFileName << "\t\thistogram: " << MCName << std::endl;
    std::cout << "  pileup Data file  :\t" << DataFileName << "\t\thistogram: " << DataName << std::endl;
    std::cout << std::endl;
    
    // PU reweighting
    LumiWeights=new edm::LumiReWeighting(MCFileName, DataFileName, MCName, DataName);
}

PileupAnalyzer::~PileupAnalyzer() {
    delete LumiWeights;
}


// ---------- PILEUP ----------

float PileupAnalyzer::GetPUWeight(const edm::Event& iEvent) {
    int nPT(0);
    // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchool2012PileupReweighting
    if(!iEvent.isRealData()) {
        edm::Handle<std::vector<PileupSummaryInfo> > PUInfo;
        iEvent.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PUInfo);
        for(std::vector<PileupSummaryInfo>::const_iterator pvi=PUInfo->begin(), pvn=PUInfo->end(); pvi!=pvn; ++pvi) {
            if(pvi->getBunchCrossing()==0) nPT=pvi->getTrueNumInteractions(); // getPU_NumInteractions();
        }
        return LumiWeights->weight( nPT );
    }
    return 1.;
}

