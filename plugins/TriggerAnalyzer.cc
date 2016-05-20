#include "TriggerAnalyzer.h"


TriggerAnalyzer::TriggerAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    TriggerToken(CColl.consumes<edm::TriggerResults>(PSet.getParameter<edm::InputTag>("trigger")))
{
    std::cout << " --- TriggerAnalyzer initialization ---" << std::endl;
    
    std::cout << std::endl;
}

TriggerAnalyzer::~TriggerAnalyzer() {

}



// ---------- TRIGGER ----------

int TriggerAnalyzer::FillTriggerBitmap(const edm::Event& iEvent, std::vector<std::string> Vect) {
    unsigned int n=Vect.size();
    std::vector<bool> Fired(n);

    edm::Handle<edm::TriggerResults> hltTriggerResults;
    iEvent.getByToken(TriggerToken, hltTriggerResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltTriggerResults);
    // Get Trigger index
    for(unsigned int i=0, in=trigNames.size(); i<in; i++) {
      for(unsigned int j=0; j<n; j++) {
        if(trigNames.triggerName(i).find(Vect[j])!=std::string::npos) {
          unsigned int index=trigNames.triggerIndex(trigNames.triggerName(i));
          Fired[j]=hltTriggerResults->accept(index);
        }
      }
    }
    
    int bitmap(0);
    for(unsigned int j=0; j<n; j++) bitmap+=pow(2, j);
    return bitmap;
}

