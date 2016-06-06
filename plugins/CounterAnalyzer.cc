#include "CounterAnalyzer.h"



CounterAnalyzer::CounterAnalyzer(const edm::ParameterSet& iConfig):
    LheToken(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheProduct")))
{
    //now do what ever initialization is needed
    usesResource("TFileService");
    Hist = fs->make<TH1F>("c_nEvents", "Event Counter", 1, 0., 1.);
    Hist->Sumw2();
    
    std::cout << " --- CounterAnalyzer initialization ---" << std::endl;
    std::cout << std::endl;
}


CounterAnalyzer::~CounterAnalyzer() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void CounterAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    float weight(1.);
    if(!iEvent.isRealData()) {
        // Declare and open collection
        edm::Handle<LHEEventProduct> LheEventCollection;
        iEvent.getByToken(LheToken, LheEventCollection);
        weight = LheEventCollection.product()->originalXWGTUP();
        weight = weight > 0. ? 1. : -1.;
    }
    Hist->Fill(0., weight);
}


// ------------ method called once each job just before starting event loop  ------------
void CounterAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void CounterAnalyzer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void CounterAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CounterAnalyzer);
