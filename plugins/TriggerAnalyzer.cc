#include "TriggerAnalyzer.h"


TriggerAnalyzer::TriggerAnalyzer() {
  
    std::cout << " - TriggerAnalyzer initialized" << std::endl;
}

TriggerAnalyzer::~TriggerAnalyzer() {

}



// ---------- TRIGGER ----------

int TriggerAnalyzer::FillTriggerBitmap(const edm::Event& iEvent, std::vector<std::string> Vect) {
    unsigned int n=Vect.size();
  //  std::vector<bool> Fired(n);

  //  edm::Handle<edm::TriggerResults> hltTriggerResults;
  //  iEvent.getByLabel(edm::InputTag("TriggerResults"), hltTriggerResults);
  //  const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltTriggerResults);
  //  // Get Trigger index
  //  for(unsigned int i=0, in=trigNames.size(); i<in; i++) {
  //    for(unsigned int j=0; j<n; j++) {
  //      if(trigNames.triggerName(i).find(Vect[j])!=std::string::npos) {
  //        unsigned int index=trigNames.triggerIndex(trigNames.triggerName(i));
  //        Fired[j]=hltTriggerResults->accept(index);
  //      }
  //    }
  //  }
  //  
    int bitmap(0);
  //  for(unsigned int j=0; j<n; j++) bitmap+=pow(2, j);
    return bitmap;
}

