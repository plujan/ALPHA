#include "PhotonAnalyzer.h"
  
PhotonAnalyzer::PhotonAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    PhotonToken(CColl.consumes<std::vector<pat::Photon> >(PSet.getParameter<edm::InputTag>("photons"))),
    VertexToken(CColl.consumes<reco::VertexCollection>(PSet.getParameter<edm::InputTag>("vertices"))),
    Photon1Id(PSet.getParameter<int>("photon1id")),
    Photon2Id(PSet.getParameter<int>("photon2id")),
    Photon1Iso(PSet.getParameter<int>("photon1iso")),
    Photon2Iso(PSet.getParameter<int>("photon2iso")),
    Photon1Pt(PSet.getParameter<double>("photon1pt")),
    Photon2Pt(PSet.getParameter<double>("photon2pt"))
{
  
    std::cout << " - PhotonAnalyzer initialized:" << std::endl;
    std::cout << "Id  :\t" << Photon1Id << "\t" << Photon2Id << std::endl;
    std::cout << "Iso :\t" << Photon1Iso << "\t" << Photon2Iso << std::endl;
    std::cout << "pT  :\t" << Photon1Pt << "\t" << Photon2Pt << std::endl;
}

PhotonAnalyzer::~PhotonAnalyzer() {

}





std::vector<pat::Photon> PhotonAnalyzer::FillPhotonVector(const edm::Event& iEvent) {
    int IdTh(Photon1Id), IsoTh(Photon1Iso);
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
    
    // Loop on Photon collection
    for(std::vector<pat::Photon>::const_iterator it=PhoCollection->begin(); it!=PhoCollection->end(); ++it) {
        if(Vect.size()>0) {
            IdTh=Photon2Id;
            IsoTh=Photon2Iso;
            PtTh=Photon2Pt;
        }
        pat::Photon ph=*it;
        // Pt and eta
        if(ph.pt()<PtTh || fabs(ph.eta())>2.5) continue;
        float pfIso = ( ph.chargedHadronIso() + std::max(ph.neutralHadronIso() + ph.photonIso() - 0.5*ph.puChargedHadronIso(), 0.) ) / ph.pt();
        ph.addUserFloat("PFIso", pfIso);
        if(IsoTh==1 && pfIso>0.15) continue;
        // Photon Id
        if(IdTh==1 && !isLoosePhoton(ph, vertex)) continue;
        Vect.push_back(ph); // Fill vector
    }
    return Vect;
}

bool PhotonAnalyzer::isLoosePhoton(pat::Photon& el, const reco::Vertex* vertex) {
    return true;
}


