// -*- C++ -*-
//
// Package:    Analysis/HZZ
// Class:      HZZ
// 
/**\class HZZ HZZ.cc Analysis/HZZ/plugins/HZZ.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alberto Zucchetta
//         Created:  Thu, 28 Apr 2016 08:28:54 GMT
//
//

#include "HZZ.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HZZ::HZZ(const edm::ParameterSet& iConfig):
    GenPSet(iConfig.getParameter<edm::ParameterSet>("genSet")),
    PileupPSet(iConfig.getParameter<edm::ParameterSet>("pileupSet")),
    TriggerPSet(iConfig.getParameter<edm::ParameterSet>("triggerSet")),
    ElectronPSet(iConfig.getParameter<edm::ParameterSet>("electronSet")),
    MuonPSet(iConfig.getParameter<edm::ParameterSet>("muonSet")),
    TauPSet(iConfig.getParameter<edm::ParameterSet>("tauSet")),
    PhotonPSet(iConfig.getParameter<edm::ParameterSet>("photonSet")),
    JetPSet(iConfig.getParameter<edm::ParameterSet>("jetSet")),
    WriteNElectrons(iConfig.getParameter<int>("writeNElectrons")),
    WriteNMuons(iConfig.getParameter<int>("writeNMuons")),
    WriteNLeptons(iConfig.getParameter<int>("writeNLeptons")),
    WriteNTaus(iConfig.getParameter<int>("writeNTaus")),
    WriteNPhotons(iConfig.getParameter<int>("writeNPhotons")),
    WriteNJets(iConfig.getParameter<int>("writeNJets")),
    HistFile(iConfig.getParameter<std::string>("histFile")),
    Verbose(iConfig.getParameter<bool>("verbose"))
{
    //now do what ever initialization is needed
    usesResource("TFileService");
    
    // Initialize Objects
    theGenAnalyzer=new GenAnalyzer(GenPSet, consumesCollector());
    thePileupAnalyzer=new PileupAnalyzer(PileupPSet, consumesCollector());
    theTriggerAnalyzer=new TriggerAnalyzer(TriggerPSet, consumesCollector());
    theElectronAnalyzer=new ElectronAnalyzer(ElectronPSet, consumesCollector());
    theMuonAnalyzer=new MuonAnalyzer(MuonPSet, consumesCollector());
    theTauAnalyzer=new TauAnalyzer(TauPSet, consumesCollector());
    thePhotonAnalyzer=new PhotonAnalyzer(PhotonPSet, consumesCollector());
    theJetAnalyzer=new JetAnalyzer(JetPSet, consumesCollector());
    //theBTagAnalyzer=new BTagAnalyzer(BTagAlgo);
    
    
    // ---------- Plots Initialization ----------
    TFileDirectory allDir=fs->mkdir("All/");
    TFileDirectory genDir=fs->mkdir("Gen/");
    TFileDirectory eleDir=fs->mkdir("Electrons/");
    TFileDirectory muoDir=fs->mkdir("Muons/");
    TFileDirectory jetDir=fs->mkdir("Jets/");
    TFileDirectory kinDir=fs->mkdir("Kin/");
    
    // Make TH1F
    std::vector<std::string> nLabels={"Trigger", "Lep #geq 2", "Z cand ", "Jets #geq 2", "Z mass ", "bJets #geq 1", "bJets #geq 2", "h mass ", "#slash{E}_{T}", "Final"};
    
    int nbins;
    float min, max;
    std::string name, title, opt;
    
    ifstream histFile(HistFile);
    if(!histFile.is_open()) {
        throw cms::Exception("HZZ Analyzer", HistFile + " file not found");
    }
    while(histFile >> name >> title >> nbins >> min >> max >> opt) {
        if(name.find('#')==std::string::npos) {
            while(title.find("~")!=std::string::npos) title=title.replace(title.find("~"), 1, " "); // Remove ~
            if(name.substr(0, 2)=="a_") Hist[name] = allDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max); //.substr(2)
            if(name.substr(0, 2)=="g_") Hist[name] = genDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            if(name.substr(0, 2)=="e_") Hist[name] = eleDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            if(name.substr(0, 2)=="m_") Hist[name] = muoDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            if(name.substr(0, 2)=="j_") Hist[name] = jetDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            if(name.substr(0, 2)=="k_") Hist[name] = kinDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            Hist[name]->Sumw2();
            Hist[name]->SetOption(opt.c_str());
            // Particular histograms
            if(name=="a_nEvents" || name=="e_nEvents" || name=="m_nEvents") for(unsigned int i=0; i<nLabels.size(); i++) Hist[name]->GetXaxis()->SetBinLabel(i+1, nLabels[i].c_str());
        }
    }
    histFile.close();
    
    std::cout << "---------- STARTING ----------" << std::endl;
}


HZZ::~HZZ() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    std::cout << "---------- ENDING  ----------" << std::endl;
    
    
    delete theGenAnalyzer;
    delete thePileupAnalyzer;
    delete theTriggerAnalyzer;
    delete theElectronAnalyzer;
    delete theMuonAnalyzer;
    delete theTauAnalyzer;
    delete thePhotonAnalyzer;
    delete theJetAnalyzer;
    //delete theBTagAnalyzer;
    
}


//
// member functions
//

// ------------ method called for each event  ------------
void HZZ::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    isMC = iEvent.isRealData();
    EventNumber = iEvent.id().event();
    LumiNumber = iEvent.luminosityBlock();
    RunNumber = iEvent.id().run();
    
    EventWeight = PUWeight = TriggerWeight = LeptonWeight = 1.;
    
    AddFourMomenta addP4;
    
    // Initialize types
    for(int i = 0; i < WriteNElectrons; i++) ObjectsFormat::ResetLeptonType(Electrons[i]);
    for(int i = 0; i < WriteNMuons; i++) ObjectsFormat::ResetLeptonType(Muons[i]);
    for(int i = 0; i < WriteNLeptons; i++) ObjectsFormat::ResetLeptonType(Leptons[i]);
    for(int i = 0; i < WriteNTaus; i++) ObjectsFormat::ResetTauType(Taus[i]);
    for(int i = 0; i < WriteNPhotons; i++) ObjectsFormat::ResetPhotonType(Photons[i]);
    for(int i = 0; i < WriteNJets; i++) ObjectsFormat::ResetJetType(Jets[i]);
    ObjectsFormat::ResetMEtType(MEt);
    ObjectsFormat::ResetCandidateType(Z2);
    ObjectsFormat::ResetCandidateType(Z1);
    
    
    // Electrons
    std::vector<pat::Electron> ElecVect = theElectronAnalyzer->FillElectronVector(iEvent);
    // Muons
    std::vector<pat::Muon> MuonVect = theMuonAnalyzer->FillMuonVector(iEvent);
    // Taus
    std::vector<pat::Tau> TauVect = theTauAnalyzer->FillTauVector(iEvent);
    // Photons
    std::vector<pat::Photon> PhotonVect = thePhotonAnalyzer->FillPhotonVector(iEvent);
    // Jets
    std::vector<pat::Jet> JetsVect = theJetAnalyzer->FillJetVector(iEvent);
    theJetAnalyzer->CleanJetsFromMuons(JetsVect, MuonVect);
    theJetAnalyzer->CleanJetsFromElectrons(JetsVect, ElecVect);
    // Missing Energy
    pat::MET MET = theJetAnalyzer->FillMetVector(iEvent);
    
    // Gen weights
    std::map<std::string, float> GenWeight = theGenAnalyzer->FillWeightsMap(iEvent);
    EventWeight *= GenWeight["event"];
    // Gen Particles
    std::vector<reco::GenParticle> GenPVect = theGenAnalyzer->FillGenVector(iEvent);
    // Gen candidates
    reco::Candidate* theGenZ = theGenAnalyzer->FindGenParticle(GenPVect, 23);
    
    // PU weight
    PUWeight = thePileupAnalyzer->GetPUWeight(iEvent);
    EventWeight *= PUWeight;
    
    // Trigger
    std::map<std::string, bool> TriggerMap = theTriggerAnalyzer->FillTriggerMap(iEvent);
    EventWeight *= TriggerWeight;
    
    
    Hist["a_nEvents"]->Fill(1., EventWeight);
    Hist["e_nEvents"]->Fill(1., EventWeight);
    Hist["m_nEvents"]->Fill(1., EventWeight);
    
    // ---------- Do analysis selections ----------
    
    // ---------- Z TO LEPTONS ----------
    bool isZtoMM(false), isZtoEE(false);
    
    if(MuonVect.size()>=2 && ElecVect.size()>=2) {
        if(MuonVect.at(0).pt() > ElecVect.at(0).pt()) isZtoMM=true;
        else isZtoEE=true;
    }
    else if(ElecVect.size()>=2) isZtoEE=true;
    else if(MuonVect.size()>=2) isZtoMM=true;
    else {if(Verbose) std::cout << " - No Iso SF OS Leptons" << std::endl; return;}
    
    // AZh lepton choice
//    int l1(0), l2(-1);
//    if(isZtoEE) {
//        for(unsigned int i=1; i<ElecVect.size(); i++) if(l2<0 && ElecVect.at(i).charge()!=ElecVect.at(l1).charge()) l2=i;
//    }
//    else if(isZtoMM) {
//        for(unsigned int i=1; i<MuonVect.size(); i++) if(l2<0 && MuonVect.at(i).charge()!=MuonVect.at(l1).charge()) l2=i;
//    }
//    if(l1<0 || l2<0) {if(Verbose) std::cout << " - No OS SF leptons" << std::endl; return;}

    // Reconstruct Z candidate
    pat::CompositeCandidate theZ2;
    if(isZtoMM) {
        theZ2.addDaughter(MuonVect.at(0));
        theZ2.addDaughter(MuonVect.at(1));
    }
    else {
        theZ2.addDaughter(ElecVect.at(0));
        theZ2.addDaughter(ElecVect.at(1));
    }
    addP4.set(theZ2);
    
    // ---------- Z TO HADRONS ----------
    pat::CompositeCandidate theZ1;
    
    if(JetsVect.size() < 2) {if(Verbose) std::cout << " - N jets < 2" << std::endl;} // return;}
    else {
        theZ1.addDaughter(JetsVect.at(0));
        theZ1.addDaughter(JetsVect.at(1));
        addP4.set(theZ1);
    }
    
    // Global candidate
    pat::CompositeCandidate theX;
    theX.addDaughter(theZ1);
    theX.addDaughter(theZ2);
    addP4.set(theX);

    // ---------- Print Summary ----------
    if(Verbose) {
        std::cout << " --- Event n. " << iEvent.id().event() << ", lumi " << iEvent.luminosityBlock() << ", run " << iEvent.id().run() << ", weight " << EventWeight << std::endl;
        std::cout << "number of electrons: " << ElecVect.size() << std::endl;
        for(unsigned int i = 0; i < ElecVect.size(); i++) std::cout << "  electron [" << i << "]\tpt: " << ElecVect[i].pt() << "\teta: " << ElecVect[i].eta() << "\tphi: " << ElecVect[i].phi() << std::endl;
        std::cout << "number of muons:     " << MuonVect.size() << std::endl;
        for(unsigned int i = 0; i < MuonVect.size(); i++) std::cout << "  muon     [" << i << "]\tpt: " << MuonVect[i].pt() << "\teta: " << MuonVect[i].eta() << "\tphi: " << MuonVect[i].phi() << std::endl;
        std::cout << "number of taus:  " << TauVect.size() << std::endl;
        for(unsigned int i = 0; i < TauVect.size(); i++) std::cout << "  tau  [" << i << "]\tpt: " << TauVect[i].pt() << "\teta: " << TauVect[i].eta() << "\tphi: " << TauVect[i].phi() << std::endl;
        std::cout << "number of photons:  " << PhotonVect.size() << std::endl;
        for(unsigned int i = 0; i < PhotonVect.size(); i++) std::cout << "  photon  [" << i << "]\tpt: " << PhotonVect[i].pt() << "\teta: " << PhotonVect[i].eta() << "\tphi: " << PhotonVect[i].phi() << std::endl;
        std::cout << "number of AK4 jets:  " << JetsVect.size() << std::endl;
        for(unsigned int i = 0; i < JetsVect.size(); i++) std::cout << "  AK4 jet  [" << i << "]\tpt: " << JetsVect[i].pt() << "\teta: " << JetsVect[i].eta() << "\tphi: " << JetsVect[i].phi() << std::endl;
        std::cout << "Missing energy:      " << MET.pt() << std::endl;
        std::cout << "Z leptonic mass:     " << theZ2.mass() << ", generated: " << (theGenZ ? theGenZ->mass() : -1.) << std::endl;
        std::cout << "Z hadronic mass:     " << theZ1.mass() << std::endl;
        std::cout << "X candidate mass:    " << theX.mass() << std::endl;
    }

    
    // ---------- Fill objects ----------
    if(isZtoEE) for(unsigned int i = 0; i < Leptons.size() && i < ElecVect.size(); i++) ObjectsFormat::FillElectronType(Leptons[i], &ElecVect[i], isMC);
    else for(unsigned int i = 0; i < Leptons.size() && i < MuonVect.size(); i++) ObjectsFormat::FillMuonType(Leptons[i], &MuonVect[i], isMC);
    for(unsigned int i = 0; i < Taus.size() && i < TauVect.size(); i++) ObjectsFormat::FillTauType(Taus[i], &TauVect[i], isMC);
    for(unsigned int i = 0; i < Photons.size() && i < PhotonVect.size(); i++) ObjectsFormat::FillPhotonType(Photons[i], &PhotonVect[i], isMC);
    for(unsigned int i = 0; i < Jets.size() && i < JetsVect.size(); i++) ObjectsFormat::FillJetType(Jets[i], &JetsVect[i], isMC);
    ObjectsFormat::FillMEtType(MEt, &MET, isMC);
    ObjectsFormat::FillCandidateType(Z2, &theZ2, isMC);
    ObjectsFormat::FillCandidateType(Z1, &theZ1, isMC);
    ObjectsFormat::FillCandidateType(X, &theX, isMC);
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
//        std::cout << "\tReconstructed Z candidate from " << (isZtoMM ? "muons" : "electrons") << " " << l1 << " and " << l2 << " with mass: " << theZ2.mass() << std::endl;
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
//    theJetAnalyzer->ApplyRecoilCorrections(MET, &MET.genMET()->p4(), &theZ2.p4(), 0);

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
void HZZ::beginJob() {
    
    // Object objects are created only one in the begin job. The reference passed to the branch has to be the same
    for(int i = 0; i < WriteNElectrons; i++) Electrons.push_back( LeptonType() );
    for(int i = 0; i < WriteNMuons; i++) Muons.push_back( LeptonType() );
    for(int i = 0; i < WriteNLeptons; i++) Leptons.push_back( LeptonType() );
    for(int i = 0; i < WriteNTaus; i++) Taus.push_back( TauType() );
    for(int i = 0; i < WriteNPhotons; i++) Photons.push_back( PhotonType() );
    for(int i = 0; i < WriteNJets; i++) Jets.push_back( JetType() );
    
    
    // Create Tree and set Branches
    tree=fs->make<TTree>("tree", "tree");
    tree->Branch("isMC", &isMC, "isMC/O");
    tree->Branch("EventNumber", &EventNumber, "EventNumber/L");
    tree->Branch("LumiNumber", &LumiNumber, "LumiNumber/L");
    tree->Branch("RunNumber", &RunNumber, "RunNumber/L");
    tree->Branch("EventWeight", &EventWeight, "EventWeight/F");
    tree->Branch("PUWeight", &PUWeight, "PUWeight/F");
    tree->Branch("TriggerWeight", &TriggerWeight, "TriggerWeight/F");
    tree->Branch("LeptonWeight", &LeptonWeight, "LeptonWeight/F");
    
    // Set Branches for objects
    for(int i = 0; i < WriteNElectrons; i++) tree->Branch(("Electron"+std::to_string(i+1)).c_str(), &(Electrons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNMuons; i++) tree->Branch(("Muon"+std::to_string(i+1)).c_str(), &(Muons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNLeptons; i++) tree->Branch(("Lepton"+std::to_string(i+1)).c_str(), &(Leptons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNTaus; i++) tree->Branch(("Tau"+std::to_string(i+1)).c_str(), &(Taus[i]), ObjectsFormat::ListTauType().c_str());
    for(int i = 0; i < WriteNPhotons; i++) tree->Branch(("Photon"+std::to_string(i+1)).c_str(), &(Photons[i]), ObjectsFormat::ListPhotonType().c_str());
    for(int i = 0; i < WriteNJets; i++) tree->Branch(("Jet"+std::to_string(i+1)).c_str(), &(Jets[i]), ObjectsFormat::ListJetType().c_str());
    tree->Branch("MEt", &MEt, ObjectsFormat::ListMEtType().c_str());
    tree->Branch("Z2", &Z2, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("Z1", &Z1, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("X", &X, ObjectsFormat::ListCandidateType().c_str());
}

// ------------ method called once each job just after ending the event loop  ------------
void HZZ::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HZZ::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HZZ);
