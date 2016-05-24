#include "TauAnalyzer.h"
  
TauAnalyzer::TauAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    TauToken(CColl.consumes<std::vector<pat::Tau> >(PSet.getParameter<edm::InputTag>("taus"))),
    VertexToken(CColl.consumes<reco::VertexCollection>(PSet.getParameter<edm::InputTag>("vertices"))),
    TauId(PSet.getParameter<int>("tauid")),
    TauPt(PSet.getParameter<double>("taupt"))

{

    std::cout << " - TauAnalyzer initialized:" << std::endl;
    std::cout << "Id  :\t" << TauId << std::endl;
    std::cout << "pT  :\t" << TauPt << std::endl;


}

TauAnalyzer::~TauAnalyzer() {
}



std::vector<pat::Tau> TauAnalyzer::FillTauVector(const edm::Event& iEvent) {
    //bool isMC(!iEvent.isRealData());
    int IdTh(TauId);
    float PtTh(TauPt);
    std::vector<pat::Tau> Vect;
    // Declare and open collection
    edm::Handle<std::vector<pat::Tau> > TauCollection;
    iEvent.getByToken(TauToken, TauCollection);
        
    edm::Handle<reco::VertexCollection> PVCollection;
    iEvent.getByToken(VertexToken, PVCollection);
    const reco::Vertex* vertex=&PVCollection->front();
    
    // Loop on Tau collection
    for(std::vector<pat::Tau>::const_iterator it=TauCollection->begin(); it!=TauCollection->end(); ++it) {
        pat::Tau tau=*it;
        // Pt and eta
        if(tau.pt()<PtTh || fabs(tau.eta())>2.5) continue;
        float pfIso = ( tau.chargedHadronIso() + std::max(tau.neutralHadronIso() + tau.photonIso() - 0.5*tau.puChargedHadronIso(), 0.) ) / tau.pt();

        //Tau CutBased and HEEP ID 2015-2016, https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2  
        //bool isPassLoose = (*LooseIdDecisions)[tauRef];
        //bool isPassMedium = (*MediumIdDecisions)[tauRef];
        //bool isPassTight = (*TightIdDecisions)[tauRef];
        //bool isPassMVANonTrigMedium = (*MVANonTrigMediumIdDecisions)[tauRef];

        //if(IdTh==1 && !isPassLoose) continue;
        //if(IdTh==2 && !isPassMedium) continue;
        //if(IdTh==3 && !isPassTight) continue;
        //if(IdTh==4 && !isPassMVANonTrigMedium) continue;

        //Fill user float
        tau.addUserFloat("pfIso", pfIso);
        //tau.addUserInt("isLoose", isPassLoose ? 1 : 0);
        //tau.addUserInt("isMedium", isPassMedium ? 1 : 0);
        //tau.addUserInt("isTight", isPassTight ? 1 : 0);        
      
        Vect.push_back(tau); // Fill vector
    }
    return Vect;
}

