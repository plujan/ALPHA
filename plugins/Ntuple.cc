// -*- C++ -*-
//
// Package:    Analysis/Ntuple
// Class:      Ntuple
// 
/**\class Ntuple Ntuple.cc Analysis/Ntuple/plugins/Ntuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alberto Zucchetta
//         Created:  Thu, 28 Apr 2016 08:28:54 GMT
//
//

#include "Ntuple.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Ntuple::Ntuple(const edm::ParameterSet& iConfig):
    ElectronPSet(iConfig.getParameter<edm::ParameterSet>("electronSet")),
    MuonPSet(iConfig.getParameter<edm::ParameterSet>("muonSet")),
    PhotonPSet(iConfig.getParameter<edm::ParameterSet>("photonSet")),
    JetPSet(iConfig.getParameter<edm::ParameterSet>("jetSet")),
    WriteNElectrons(iConfig.getParameter<int>("writeNElectrons")),
    WriteNMuons(iConfig.getParameter<int>("writeNMuons")),
    WriteNLeptons(iConfig.getParameter<int>("writeNLeptons")),
    WriteNJets(iConfig.getParameter<int>("writeNJets")),
    Verbose(iConfig.getParameter<bool>("verbose"))
{
    //now do what ever initialization is needed
    usesResource("TFileService");
    
    // Initialize Objects
    theGenAnalyzer=new GenAnalyzer();
    thePileupAnalyzer=new PileupAnalyzer();
    theTriggerAnalyzer=new TriggerAnalyzer();
    theElectronAnalyzer=new ElectronAnalyzer(ElectronPSet, consumesCollector());
    theMuonAnalyzer=new MuonAnalyzer(MuonPSet, consumesCollector());
    thePhotonAnalyzer=new PhotonAnalyzer(PhotonPSet, consumesCollector());
    theJetAnalyzer=new JetAnalyzer(JetPSet, consumesCollector());
    //theBTagAnalyzer=new BTagAnalyzer(BTagAlgo);
    
    std::cout << "---------- STARTING ----------" << std::endl;
}


Ntuple::~Ntuple() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    std::cout << "---------- ENDING  ----------" << std::endl;
    
    
    delete theGenAnalyzer;
    delete thePileupAnalyzer;
    delete theTriggerAnalyzer;
    delete theElectronAnalyzer;
    delete theMuonAnalyzer;
    delete thePhotonAnalyzer;
    delete theJetAnalyzer;
    //delete theBTagAnalyzer;
    
}


//
// member functions
//

// ------------ method called for each event  ------------
void Ntuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    isMC = iEvent.isRealData();
    EventNumber = iEvent.id().event();
    LumiNumber = iEvent.luminosityBlock();
    RunNumber = iEvent.id().run();
    
    float EventWeight(1.), PUWeight(1.), TriggerWeight(1.), LeptonWeight(1.);
    
    // Initialize types
    for(int i = 0; i < WriteNElectrons; i++) ObjectsFormat::ResetLeptonType(Electrons[i]);
    for(int i = 0; i < WriteNMuons; i++) ObjectsFormat::ResetLeptonType(Muons[i]);
    for(int i = 0; i < WriteNLeptons; i++) ObjectsFormat::ResetLeptonType(Leptons[i]);
    for(int i = 0; i < WriteNJets; i++) ObjectsFormat::ResetJetType(Jets[i]);
    
    // Electrons
    std::vector<pat::Electron> ElecVect=theElectronAnalyzer->FillElectronVector(iEvent);
    // Muons
    std::vector<pat::Muon> MuonVect=theMuonAnalyzer->FillMuonVector(iEvent);
    // Jets
    std::vector<pat::Jet> JetsVect=theJetAnalyzer->FillJetVector(iEvent);
    // Missing Energy
    pat::MET MEt = theJetAnalyzer->FillMetVector(iEvent);
    
    // PU weight
    PUWeight=theGenAnalyzer->GetPUWeight(iEvent);
    EventWeight*=PUWeight;
    
    // Trigger
    std::vector<std::string> TrigNames;
    TrigNames.push_back("HLT_Mu17_Mu8");
    TrigNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL");
    int TrigBit=theTriggerAnalyzer->FillTriggerBitmap(iEvent, TrigNames);
    
    
    // ---------- Print Summary ----------
    if(Verbose) {
        std::cout << " --- Event n. " << iEvent.id().event() << ", lumi " << iEvent.luminosityBlock() << ", run " << iEvent.id().run() << std::endl;
        std::cout << "number of electrons: " << ElecVect.size() << std::endl;
        for(unsigned int i = 0; i < ElecVect.size(); i++) std::cout << "  electron [" << i << "]\tpt: " << ElecVect[i].pt() << "\teta: " << ElecVect[i].eta() << "\tphi: " << ElecVect[i].phi() << std::endl;
        std::cout << "number of muons:     " << MuonVect.size() << std::endl;
        for(unsigned int i = 0; i < MuonVect.size(); i++) std::cout << "  muon     [" << i << "]\tpt: " << MuonVect[i].pt() << "\teta: " << MuonVect[i].eta() << "\tphi: " << MuonVect[i].phi() << std::endl;
        std::cout << "number of AK4 jets:  " << JetsVect.size() << std::endl;
        for(unsigned int i = 0; i < JetsVect.size(); i++) std::cout << "  AK4 jet  [" << i << "]\tpt: " << JetsVect[i].pt() << "\teta: " << JetsVect[i].eta() << "\tphi: " << JetsVect[i].phi() << std::endl;
        
        if(ElecVect.size() > MuonVect.size()) {
            for(unsigned int i = 0; i < Leptons.size() && i < ElecVect.size(); i++) ObjectsFormat::FillElectronType(Leptons[i], &ElecVect[i], isMC);
        }
        else {
            for(unsigned int i = 0; i < Leptons.size() && i < MuonVect.size(); i++) ObjectsFormat::FillMuonType(Leptons[i], &MuonVect[i], isMC);
        }
        for(unsigned int i = 0; i < Jets.size() && i < JetsVect.size(); i++) ObjectsFormat::FillJetType(Jets[i], &JetsVect[i], isMC);
    }
    
    // ---------- Z TO LEPTONS ----------
//    bool isZtoMM(false), isZtoEE(false);
//    int l1(0), l2(-1);
//    
//    if(MuonVect.size()>=2 && ElecVect.size()>=2) {
//        if(MuonVect.at(0).pt() > ElecVect.at(0).pt()) {isZtoMM=true; isZtoEE=false;}
//        else {isZtoMM=false; isZtoEE=true;}
//    }
//    else if(ElecVect.size()>=2) {isZtoMM=false; isZtoEE=true;}
//    else if(MuonVect.size()>=2) {isZtoMM=true; isZtoEE=false;}
//    else {if(Verbose) std::cout << " - No Iso SF OS Leptons" << std::endl; return;}

//    if(isZtoEE) {
//        for(unsigned int i=1; i<ElecVect.size(); i++) if(l2<0 && ElecVect.at(i).charge()!=ElecVect.at(l1).charge()) l2=i;
//    }
//    else {
//        for(unsigned int i=1; i<MuonVect.size(); i++) if(l2<0 && MuonVect.at(i).charge()!=MuonVect.at(l1).charge()) l2=i;
//    }
//    if(l1<0 || l2<0) {if(Verbose) std::cout << " - No OS SF leptons" << std::endl; return;}
//    
//    if(Verbose) {
//        std::cout << "\tEvent n. " << iEvent.id().event() << ", weight: " << EventWeight << std::endl;
//        std::cout << "\tTrigger fired: " << TrigBit << std::endl;
//        //for(unsigned int i=0; i<TrigNames.size(); i++) if(TrigBit%==0) std::cout << "\t\t" << i << " - " << TrigNames[i] << std::endl;
//        std::cout << "\tNumber of Electrons: " << ElecVect.size() << std::endl;
//        for(unsigned int i=0; i<ElecVect.size(); i++) std::cout << "\t\t" << i << " - pt: " << ElecVect[i].pt() << ", eta: " << ElecVect[i].eta() << ", phi: " << ElecVect[i].phi() << std::endl;
//        std::cout << "\tNumber of Muons: " << MuonVect.size() << std::endl;
//        for(unsigned int i=0; i<MuonVect.size(); i++) std::cout << "\t\t" << i << " - pt: " << MuonVect[i].pt() << ", eta: " << MuonVect[i].eta() << ", phi: " << MuonVect[i].phi() << std::endl;
//        std::cout << "\tNumber of Jets: " << JetsVect.size() << std::endl;
//        for(unsigned int i=0; i<JetsVect.size(); i++) std::cout << "\t\t" << i << " - pt: " << JetsVect[i].pt() << ", eta: " << JetsVect[i].eta() << ", phi: " << JetsVect[i].phi() << std::endl;
//        std::cout << "\tMissing Energy: " << MEt.pt() << std::endl;
//        std::cout << std::endl;
//    }
//    
//    // Reconstruct Z candidate
//    pat::CompositeCandidate theZ;
//    if(isZtoEE) {
//        theZ.addDaughter(ElecVect.at(l1));
//        theZ.addDaughter(ElecVect.at(l2));
//    }
//    else {
//        theZ.addDaughter(MuonVect.at(l1));
//        theZ.addDaughter(MuonVect.at(l2));
//    }
//    AddFourMomenta addP4;
//    addP4.set(theZ);
//    if(theZ.mass()<50.) {if(Verbose) std::cout << " - Z off-shell" << std::endl; return;}
//    
//    // Lepton and Trigger SF
//    if(isMC) {
//        if(isZtoEE) {
//            TriggerWeight*=theElectronAnalyzer->GetDoubleElectronTriggerSF(ElecVect.at(l1), ElecVect.at(l2));
//            LeptonWeight*=theElectronAnalyzer->GetElectronIdSF(ElecVect.at(l1));
//            LeptonWeight*=theElectronAnalyzer->GetElectronIdSF(ElecVect.at(l2));
//            LeptonWeight*=theElectronAnalyzer->GetElectronIsoSF(ElecVect.at(l1));
//            LeptonWeight*=theElectronAnalyzer->GetElectronIsoSF(ElecVect.at(l2));
//        }
//        else {
//            TriggerWeight*=theMuonAnalyzer->GetDoubleMuonTriggerSF(MuonVect.at(l1), MuonVect.at(l2));
//            LeptonWeight*=theMuonAnalyzer->GetMuonIdSF(MuonVect.at(l1));
//            LeptonWeight*=theMuonAnalyzer->GetMuonIdSF(MuonVect.at(l2));
//            LeptonWeight*=theMuonAnalyzer->GetMuonIsoSF(MuonVect.at(l1));
//            LeptonWeight*=theMuonAnalyzer->GetMuonIsoSF(MuonVect.at(l2));
//        }
//    }
//    EventWeight*=TriggerWeight;
//    EventWeight*=LeptonWeight;
//    
//    if(Verbose) {
//        std::cout << "\tReconstructed Z candidate from " << (isZtoMM ? "muons" : "electrons") << " " << l1 << " and " << l2 << " with mass: " << theZ.mass() << std::endl;
//    }
//    
//    // FatJet
//    if(JetsVect.size()<=0) {if(Verbose) std::cout << " - No Fat Jet" << std::endl; return;}
//    const pat::Jet* fatJet=&JetsVect.at(0);
//    int nSubJets=fatJet->numberOfDaughters();
//    const pat::Jet* subJet1 = dynamic_cast<const pat::Jet*>(fatJet->daughter(0));
//    const pat::Jet* subJet2 = dynamic_cast<const pat::Jet*>(fatJet->daughter(1));
//    //double subjet0Bdisc = subjet->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
//    
//    
//    if(Verbose) {
//        std::cout << "\tReconstructed " << JetsVect.size() << " fatjets:" << std::endl;
//        for(int i=0; i<(int)JetsVect.size(); i++) std::cout << "\t  fatjet " << i << " has pT: " << JetsVect[i].pt() << ", mass: " << JetsVect[i].mass() << " and CSVIVFV2: " << JetsVect[i].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
//        std::cout << "\tLeading fat jet with " << nSubJets << " subjets, leading has pT: " << fatJet->pt() << ", mass: " << fatJet->mass() << " and CSVIVFV2: " << fatJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
//        if(subJet1) std::cout << "\t  subjet 1 has pT: " << subJet1->pt() << ", mass: " << subJet1->mass() << " and CSVIVFV2: " << subJet1->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
//        if(subJet2) std::cout << "\t  subjet 2 has pT: " << subJet2->pt() << ", mass: " << subJet2->mass() << " and CSVIVFV2: " << subJet2->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
//    }

    tree->Fill();

}

//#ifdef THIS_IS_AN_EVENT_EXAMPLE
//   Handle<ExampleData> pIn;
//   iEvent.getByLabel("example",pIn);
//#endif
//   
//#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//   ESHandle<SetupData> pSetup;
//   iSetup.get<SetupRecord>().get(pSetup);
//#endif



// ------------ method called once each job just before starting event loop  ------------
void Ntuple::beginJob() {
    
    // Object objects are created only one in the begin job. The reference passed to the branch has to be the same
    for(int i = 0; i < WriteNElectrons; i++) Electrons.push_back( LeptonType() );
    for(int i = 0; i < WriteNMuons; i++) Muons.push_back( LeptonType() );
    for(int i = 0; i < WriteNLeptons; i++) Leptons.push_back( LeptonType() );
    for(int i = 0; i < WriteNJets; i++) Jets.push_back( JetType() );
    
    // Create Tree and set Branches
    tree=fs->make<TTree>("tree", "tree");
    tree->Branch("isMC", &isMC, "isMC/O");
    tree->Branch("EventNumber", &EventNumber, "EventNumber/L");
    tree->Branch("LumiNumber", &LumiNumber, "LumiNumber/L");
    tree->Branch("RunNumber", &RunNumber, "RunNumber/L");
    
    // Set Branches for objects
    for(int i = 0; i < WriteNElectrons; i++) tree->Branch(("Electron"+std::to_string(i+1)).c_str(), &(Electrons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNMuons; i++) tree->Branch(("Muon"+std::to_string(i+1)).c_str(), &(Muons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNLeptons; i++) tree->Branch(("Lepton"+std::to_string(i+1)).c_str(), &(Leptons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNJets; i++) tree->Branch(("Jet"+std::to_string(i+1)).c_str(), &(Jets[i]), ObjectsFormat::ListJetType().c_str());
    
}

// ------------ method called once each job just after ending the event loop  ------------
void Ntuple::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Ntuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuple);
