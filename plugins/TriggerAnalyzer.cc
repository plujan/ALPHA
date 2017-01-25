#include "TriggerAnalyzer.h"


TriggerAnalyzer::TriggerAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    TriggerToken(CColl.consumes<edm::TriggerResults>(PSet.getParameter<edm::InputTag>("trigger"))),
    TriggerList(PSet.getParameter<std::vector<std::string> >("paths")),
    MetFiltersToken(CColl.consumes<edm::TriggerResults>(PSet.getParameter<edm::InputTag>("metfilters"))),
    MetFiltersList(PSet.getParameter<std::vector<std::string> >("metpaths")),
    BadPFMuonFilterToken(CColl.consumes<bool>(PSet.getParameter<edm::InputTag>("badPFMuonFilter"))),
    BadChCandFilterToken(CColl.consumes<bool>(PSet.getParameter<edm::InputTag>("badChCandFilter")))
{
    std::cout << " --- TriggerAnalyzer initialization ---" << std::endl;
    std::cout << "  HLT paths:" << std::endl;
    for(unsigned int i = 0; i < TriggerList.size(); i++) std::cout << "    " << TriggerList[i] << "*" << std::endl;
    std::cout << "  Met filters:" << std::endl;
    for(unsigned int i = 0; i < MetFiltersList.size(); i++) std::cout << "    " << MetFiltersList[i] << "*" << std::endl;
    std::cout << std::endl;
}

TriggerAnalyzer::~TriggerAnalyzer() {

}



// ---------- TRIGGER ----------

void TriggerAnalyzer::FillTriggerMap(const edm::Event& iEvent, std::map<std::string, bool>& Map) {

    edm::Handle<edm::TriggerResults> hltTriggerResults;
    iEvent.getByToken(TriggerToken, hltTriggerResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltTriggerResults);
    
    //for(unsigned int j=0, in=trigNames.size(); j < in; j++) std::cout << trigNames.triggerName(j) << std::endl;
    
    // Get Trigger index
    for(unsigned int i = 0; i < TriggerList.size(); i++) {
        Map[TriggerList[i]] = false;
        for(unsigned int j=0, in=trigNames.size(); j < in; j++) {
            if(trigNames.triggerName(j).find(TriggerList[i]) != std::string::npos) {
                unsigned int index = trigNames.triggerIndex(trigNames.triggerName(j));
                if(hltTriggerResults->accept(index)) Map[TriggerList[i]] = true;
            }
        }
    }
}


void TriggerAnalyzer::FillMetFiltersMap(const edm::Event& iEvent, std::map<std::string, bool>& Map) { //, edm::EDGetTokenT<edm::TriggerResults>& iToken, std::vector<std::string>& iList) {

    edm::Handle<edm::TriggerResults> hltTriggerResults;
    iEvent.getByToken(MetFiltersToken, hltTriggerResults);//(iToken, hltTriggerResults);//
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltTriggerResults);
    
    //for(unsigned int j=0, in=trigNames.size(); j < in; j++) std::cout << trigNames.triggerName(j) << std::endl;
    
    // Get Trigger index
    for(unsigned int i = 0; i < MetFiltersList.size(); i++) {
        Map[MetFiltersList[i]] = false;
        for(unsigned int j=0, in=trigNames.size(); j < in; j++) {
            if(trigNames.triggerName(j).find(MetFiltersList[i]) != std::string::npos) {
                unsigned int index = trigNames.triggerIndex(trigNames.triggerName(j));
                if(hltTriggerResults->accept(index)) Map[MetFiltersList[i]] = true;
            }
        }
    }
}

bool TriggerAnalyzer::GetBadPFMuonFlag(const edm::Event& iEvent) { 

    edm::Handle<bool> ifilterbadPFMuon;
    iEvent.getByToken(BadPFMuonFilterToken, ifilterbadPFMuon);
    bool filterbadPFMuon = *ifilterbadPFMuon;

    return filterbadPFMuon;
}

bool TriggerAnalyzer::GetBadChCandFlag(const edm::Event& iEvent) { 

    edm::Handle<bool> ifilterbadChCand;
    iEvent.getByToken(BadChCandFilterToken, ifilterbadChCand);
    bool filterbadChCand = *ifilterbadChCand;

    return filterbadChCand;
}
