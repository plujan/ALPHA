#include "TriggerAnalyzer.h"


TriggerAnalyzer::TriggerAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    TriggerToken(CColl.consumes<edm::TriggerResults>(PSet.getParameter<edm::InputTag>("trigger"))),
    TriggerList(PSet.getParameter<std::vector<std::string> >("paths"))
{
    std::cout << " --- TriggerAnalyzer initialization ---" << std::endl;
    std::cout << "  HLT paths:" << std::endl;
    for(unsigned int i = 0; i < TriggerList.size(); i++) std::cout << "    " << TriggerList[i] << "*" << std::endl;
    std::cout << std::endl;
}

TriggerAnalyzer::~TriggerAnalyzer() {

}



// ---------- TRIGGER ----------

std::map<std::string, bool> TriggerAnalyzer::FillTriggerMap(const edm::Event& iEvent) {
    std::map<std::string, bool> Map;
    unsigned int n = TriggerList.size();
    std::vector<bool> Fired(n);

    edm::Handle<edm::TriggerResults> hltTriggerResults;
    iEvent.getByToken(TriggerToken, hltTriggerResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltTriggerResults);
    // Get Trigger index
    for(unsigned int i = 0; i < n; i++) {
        Map[TriggerList[i]] = false;
        for(unsigned int j=0, in=trigNames.size(); j < in; j++) {
            if(trigNames.triggerName(j).find(TriggerList[i]) != std::string::npos) {
                unsigned int index = trigNames.triggerIndex(trigNames.triggerName(j));
                if(hltTriggerResults->accept(index)) Map[TriggerList[i]] = true;
            }
        }
    }
    return Map;
}

