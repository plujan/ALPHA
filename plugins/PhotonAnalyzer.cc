#include "PhotonAnalyzer.h"
  
PhotonAnalyzer::PhotonAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    PhotonToken(CColl.consumes<std::vector<pat::Photon> >(PSet.getParameter<edm::InputTag>("photons"))),
    VertexToken(CColl.consumes<reco::VertexCollection>(PSet.getParameter<edm::InputTag>("vertices"))),
    //PhoLooseIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("phoLooseIdMap"))),
    //PhoMediumIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("phoMediumIdMap"))),
    //PhoTightIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("phoTightIdMap"))),
    //PhoMVANonTrigMediumIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("phoMVANonTrigMediumIdMap"))),
    Photon1Id(PSet.getParameter<int>("photon1id")),
    Photon2Id(PSet.getParameter<int>("photon2id")),
    //Photon1Iso(PSet.getParameter<int>("photon1iso")),
    //Photon2Iso(PSet.getParameter<int>("photon2iso")),
    Photon1Pt(PSet.getParameter<double>("photon1pt")),
    Photon2Pt(PSet.getParameter<double>("photon2pt"))
{
  
    std::cout << " - PhotonAnalyzer initialized:" << std::endl;
    std::cout << "Id  :\t" << Photon1Id << "\t" << Photon2Id << std::endl;
    //std::cout << "Iso :\t" << Photon1Iso << "\t" << Photon2Iso << std::endl;
    std::cout << "pT  :\t" << Photon1Pt << "\t" << Photon2Pt << std::endl;
}

PhotonAnalyzer::~PhotonAnalyzer() {

}





std::vector<pat::Photon> PhotonAnalyzer::FillPhotonVector(const edm::Event& iEvent) {
    bool isMC(!iEvent.isRealData());
    int IdTh(Photon1Id);//, IsoTh(Photon1Iso);
    float PtTh(Photon1Pt);
    std::vector<pat::Photon> Vect;
    // Declare and open collection
    edm::Handle<std::vector<pat::Photon> > PhoCollection;
    iEvent.getByToken(PhotonToken, PhoCollection);
    
    //edm::Handle<std::vector<pat::Conversion> > EleConv;
    //iEvent.getByToken(edm::InputTag("patConversions"), EleConv);
    
    edm::Handle<reco::VertexCollection> PVCollection;
    iEvent.getByToken(VertexToken, PVCollection);
    const reco::Vertex* vertex=&PVCollection->front();
    

    //value map for ID 2015-2016: not yet working on CMSSW 8XX
    //edm::Handle<edm::ValueMap<bool> > LooseIdDecisions;
    //edm::Handle<edm::ValueMap<bool> > MediumIdDecisions;
    //edm::Handle<edm::ValueMap<bool> > TightIdDecisions;
    //edm::Handle<edm::ValueMap<bool> > MVANonTrigMediumIdDecisions;
    //iEvent.getByToken(PhoLooseIdMapToken, LooseIdDecisions);
    //iEvent.getByToken(PhoMediumIdMapToken, MediumIdDecisions);
    //iEvent.getByToken(PhoTightIdMapToken, TightIdDecisions);
    //iEvent.getByToken(PhoMVANonTrigMediumIdMapToken, MVANonTrigMediumIdDecisions);
    unsigned int phIdx = 0;

    // Loop on Photon collection
    for(std::vector<pat::Photon>::const_iterator it=PhoCollection->begin(); it!=PhoCollection->end(); ++it) {
        if(Vect.size()>0) {
            IdTh=Photon2Id;
            //IsoTh=Photon2Iso;
            PtTh=Photon2Pt;
        }
        pat::Photon ph=*it;
	pat::PhotonRef phRef(PhoCollection,phIdx);
        // Pt and eta
        if(ph.pt()<PtTh || fabs(ph.eta())>2.5) continue;
        float pfIso = ( ph.chargedHadronIso() + std::max(ph.neutralHadronIso() + ph.photonIso() - 0.5*ph.puChargedHadronIso(), 0.) ) / ph.pt();

        //Photon CutBased and HEEP ID 2015-2016, https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2  
        //bool isPassLoose = (*LooseIdDecisions)[phRef];
        //bool isPassMedium = (*MediumIdDecisions)[phRef];
        //bool isPassTight = (*TightIdDecisions)[phRef];
        //bool isPassMVANonTrigMedium = (*MVANonTrigMediumIdDecisions)[phRef];
        //Dummy test:
        bool isDummy = true;
        bool isPassLoose = false;//(*LooseIdDecisions)[phRef];
	bool isPassMedium = false;//(*MediumIdDecisions)[phRef];
	bool isPassTight = false;//(*TightIdDecisions)[phRef];        
	bool isPassMVANonTrigMedium = false;
        if(IdTh==0 && !isDummy) continue;
        if(IdTh==1 && !isPassLoose) continue;
        if(IdTh==2 && !isPassMedium) continue;
        if(IdTh==3 && !isPassTight) continue;
        if(IdTh==4 && !isPassMVANonTrigMedium) continue;

        //Fill user float
        ph.addUserFloat("PFIso", pfIso);
        ph.addUserInt("isLoose", isPassLoose ? 1 : 0);
        ph.addUserInt("isMedium", isPassMedium ? 1 : 0);
        ph.addUserInt("isTight", isPassTight ? 1 : 0);
        ph.addUserInt("isMVANonTrigMedium", isPassMVANonTrigMedium ? 1 : 0);
        
        ++phIdx;
      
        Vect.push_back(ph); // Fill vector
    }
    return Vect;
}

/*bool PhotonAnalyzer::isLoosePhoton(pat::Photon& el, const reco::Vertex* vertex) {
    return true;
}*/


