// -*- C++ -*-
//
// Package:    Analysis/SSleptons
// Class:      SSleptons
// 
/**\class SSleptons SSleptons.cc Analysis/SSleptons/plugins/SSleptons.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alberto Zucchetta
//         Created:  Thu, 28 Apr 2016 08:28:54 GMT
//
//

#include "SSleptons.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SSleptons::SSleptons(const edm::ParameterSet& iConfig):

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
    
    std::vector<std::string> TriggerList(TriggerPSet.getParameter<std::vector<std::string> >("paths"));
    for(unsigned int i = 0; i < TriggerList.size(); i++) TriggerMap[ TriggerList[i] ] = false;
    
    std::cout << "---------- STARTING ----------" << std::endl;
}


SSleptons::~SSleptons() {
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
void SSleptons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    isMC = !iEvent.isRealData();
    EventNumber = iEvent.id().event();
    LumiNumber = iEvent.luminosityBlock();
    RunNumber = iEvent.id().run();

    //float EventWeight(1.), PUWeight(1.), TriggerWeight(1.), LeptonWeight(1.);
    EventWeight = PUWeight = TriggerWeight = LeptonWeight = 1.;
    ZewkWeight  = WewkWeight = 1.;
    nPV = 1;
    
    // Initialize types
    for(int i = 0; i < WriteNElectrons; i++) ObjectsFormat::ResetLeptonType(Electrons[i]);
    for(int i = 0; i < WriteNMuons; i++) ObjectsFormat::ResetLeptonType(Muons[i]);
    for(int i = 0; i < WriteNLeptons; i++) ObjectsFormat::ResetLeptonType(Leptons[i]);
    for(int i = 0; i < WriteNTaus; i++) ObjectsFormat::ResetTauType(Taus[i]);
    for(int i = 0; i < WriteNPhotons; i++) ObjectsFormat::ResetPhotonType(Photons[i]);
    for(int i = 0; i < WriteNJets; i++) ObjectsFormat::ResetJetType(Jets[i]);
    ObjectsFormat::ResetMEtType(MEt);
 
     // -----------------------------------
    //           READ OBJECTS
    // -----------------------------------

    
           
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
    // Missing Energy
    pat::MET MET = theJetAnalyzer->FillMetVector(iEvent);
    
    // Gen weights
    std::map<int, float> GenWeight = theGenAnalyzer->FillWeightsMap(iEvent);
    EventWeight *= GenWeight[-1];
    // Gen Particles
    std::vector<reco::GenParticle> GenPVect = theGenAnalyzer->FillGenVector(iEvent);
    // Gen candidates
    reco::Candidate* theGenZ = theGenAnalyzer->FindGenParticle(GenPVect, 23);
    reco::Candidate* theGenW = theGenAnalyzer->FindGenParticle(GenPVect, 24);
    // EWK corrections
    if(theGenZ) ZewkWeight = theGenAnalyzer->GetZewkWeight(theGenZ->pt());
    if(theGenW) WewkWeight = theGenAnalyzer->GetWewkWeight(theGenW->pt());
    EventWeight *= ZewkWeight * WewkWeight;

// ================================================
// MCgen  info 
// ================================================
    std::vector<const reco::Candidate *> GenMuonsVect;
    GenMuonsIDVect.clear(); 
    GenMuonsPtVect.clear(); 
    GenMuonsEtaVect.clear(); 
    GenMuonsPhiVect.clear(); 
    int nMuoncount = 0;
    for( unsigned int a=0; a<GenPVect.size(); a++){  //loop on GenPVect
      //std::cout << "entro nel loop di GenPVect" << std::endl;
      int istore = 0; // flag to store the particle
    // require particle to be W/Z, with correct status (= 62 ) and daughters....        
//      if( GenPVect[a].numberOfDaughters()>0 && 
//         (abs(GenPVect[a].pdgId())==24 || GenPVect[a].pdgId()==23) && GenPVect[a].status()==62) istore = 1;  // for SLtop: status = 52
      if( GenPVect[a].numberOfDaughters() > 0 ) {
         if ( abs(GenPVect[a].pdgId())==15  )                                       istore = 1; 
         if ( abs(GenPVect[a].pdgId())==24 || GenPVect[a].pdgId()==23 )             istore = 1; 
         if ( abs(GenPVect[a].pdgId()) > 200  && abs(GenPVect[a].pdgId()) < 600   ) istore = 1; 	 
	 if ( abs(GenPVect[a].pdgId()) > 4000 && abs(GenPVect[a].pdgId()) < 6000  ) istore = 1; }
        if( istore == 1) { 
	    //int ndaughters = GenPVect[a].numberOfDaughters();//dimensione generica delle figlie della W/Z
          //  std::cout << " found particle id = "<< GenPVect[a].pdgId() << " status = " << GenPVect[a].status() << std::endl;
            for(unsigned int it=0; it<GenPVect[a].numberOfDaughters(); it++) { //loop on daughters        
	      //std::cout << "ho " << GenPVect[a].numberOfDaughters() << " figlie" << std::endl;
                const reco::Candidate * madre = dynamic_cast<reco::Candidate*>(&GenPVect[a]);
                //std::cout << "fatto il dynamic cast" << std::endl;
                //for( size_t il = 0; il < madre->numberOfDaughters(); ++ il ) {
                //    std::cout << "loop sulle figlie: figlia n " << il << std::endl;
                //    const reco::Candidate * daughter = madre->daughter( il );
                //    std::cout << "daughter pt: " << daughter->pt() << std::endl;
                //}
                //std::cout << "Le figlie sono" << madre->daughter(it)->pdgId() << std::endl;
	        // save muons in GenMuonsVect
                if( abs(madre->daughter(it)->pdgId())==13){
		  if( nMuoncount < 5 ) {
                 //   std::cout << "ho delle figlie muoni" << std::endl;
		    GenMuonsIDVect[nMuoncount] = GenPVect[a].pdgId();  // store the mother...
                    GenMuonsVect.push_back( madre->daughter(it) );
		    GenMuonsPtVect.push_back( madre->daughter(it)->pt() );
		    GenMuonsEtaVect.push_back( madre->daughter(it)->eta() );
		    GenMuonsPhiVect.push_back( madre->daughter(it)->phi() );
                    nMuoncount++;
	//	    if ( abs(GenPVect[a].pdgId()) > 200  && abs(GenPVect[a].pdgId()) < 600   ) 
        //              std::cout << "gen mu pt: " << madre->daughter(it)->pt() << " mothID " <<  GenPVect[a].pdgId() << std::endl;
	//	    if ( abs(GenPVect[a].pdgId()) == 15   ) 
        //               std::cout << "gen mu pt: " << madre->daughter(it)->pt() << " mothID " <<  GenPVect[a].pdgId() << std::endl;
	         }
                } // endif muon
            }// end loop on  daughters
        }// endif  store particle 
    }// next particle in GenPVect
    nGenMuons = GenMuonsVect.size() ; // number of true muons from W/Z
 //   if( nGenMuons > 0 )  std::cout << "N = " << nGenMuons << " gen muons" << std::endl;

 
    
    // PU weight
    PUWeight = thePileupAnalyzer->GetPUWeight(iEvent);
    nPV = thePileupAnalyzer->GetPV(iEvent);
    EventWeight *= PUWeight;
    
    // Trigger
    theTriggerAnalyzer->FillTriggerMap(iEvent, TriggerMap);
    EventWeight *= TriggerWeight;
    
    
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
    }
    
    // ---------- Do analysis selections ----------
    // ...
    
    // ---------- Fill objects ----------
    if(Verbose) std::cout << " - Filling objects ----------- " << std::endl;
    // electron & muons (separated)
for(unsigned int i = 0; i < ElecVect.size(); i++) ObjectsFormat::FillElectronType(Electrons[i], &ElecVect[i], isMC);
for(unsigned int i = 0; i < MuonVect.size(); i++) ObjectsFormat::FillMuonType(Muons[i], &MuonVect[i], isMC);    
// leptons (together)
   for(unsigned int i = 0; i < Leptons.size() && i < MuonVect.size(); i++) ObjectsFormat::FillMuonType(Leptons[i], &MuonVect[i], isMC);
//   for(unsigned int i = 0;  i < MuonVect.size(); i++) ObjectsFormat::FillMuonType(Leptons[i], &MuonVect[i], isMC);

   for(unsigned int i = MuonVect.size(); i < Leptons.size() && i < ElecVect.size(); i++) ObjectsFormat::FillElectronType(Leptons[i], &ElecVect[i], isMC);

//    if(ElecVect.size() > MuonVect.size()) {
//        for(unsigned int i = 0; i < Leptons.size() && i < ElecVect.size(); i++) ObjectsFormat::FillElectronType(Leptons[i], &ElecVect[i], isMC);
//    }
//    else {
//        for(unsigned int i = 0; i < Leptons.size() && i < MuonVect.size(); i++) ObjectsFormat::FillMuonType(Leptons[i], &MuonVect[i], isMC);
//    }
    for(unsigned int i = 0; i < Taus.size() && i < TauVect.size(); i++) ObjectsFormat::FillTauType(Taus[i], &TauVect[i], isMC);
    for(unsigned int i = 0; i < Photons.size() && i < PhotonVect.size(); i++) ObjectsFormat::FillPhotonType(Photons[i], &PhotonVect[i], isMC);
    for(unsigned int i = 0; i < Jets.size() && i < JetsVect.size(); i++) ObjectsFormat::FillJetType(Jets[i], &JetsVect[i], isMC);
    ObjectsFormat::FillMEtType(MEt, &MET, isMC);
    

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
//    EventWeight *= TriggerWeight;
    EventWeight *= LeptonWeight;
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
//    theJetAnalyzer->ApplyRecoilCorrections(MET, &MET.genMET()->p4(), &theZ.p4(), 0);

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
void SSleptons::beginJob() {
    
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
    
//  cloned from Dibosons.cc =================
    tree->Branch("EventWeight", &EventWeight, "EventWeight/F");
    tree->Branch("PUWeight", &PUWeight, "PUWeight/F");
    tree->Branch("TriggerWeight", &TriggerWeight, "TriggerWeight/F");
    tree->Branch("LeptonWeight", &LeptonWeight, "LeptonWeight/F");

    // Set trigger branches
    for(auto it = TriggerMap.begin(); it != TriggerMap.end(); it++) tree->Branch(it->first.c_str(), &(it->second), (it->first+"/O").c_str());


// ==============================
   
     tree->Branch("nPV", &nPV, "nPV/I");
     
        
    // Set Branches for objects  
    
    
    // MCgen muons =========================================
    tree->Branch("nGenMuons", &nGenMuons, "nGenMuons/I");
    std::cout << "nGenMuons in begin job " << nGenMuons << std::endl;
    std::cout << "GenMuonsPtVec size in begin job " << GenMuonsPtVect.size() << std::endl;
    for(unsigned int i = 0; i < GenMuonsPtVect.size(); i++) {
       tree->Branch(("GenMuon"+std::to_string(i+1)+"ID").c_str() , &(GenMuonsIDVect[i]) , ("GenMuon"+std::to_string(i+1)+"ID/I").c_str());
       tree->Branch(("GenMuon"+std::to_string(i+1)+"pt").c_str() , &(GenMuonsPtVect[i]) , ("GenMuon"+std::to_string(i+1)+"pt/F").c_str());
       tree->Branch(("GenMuon"+std::to_string(i+1)+"eta").c_str(), &(GenMuonsEtaVect[i]), ("GenMuon"+std::to_string(i+1)+"eta/F").c_str());
       tree->Branch(("GenMuon"+std::to_string(i+1)+"phi").c_str(), &(GenMuonsPhiVect[i]), ("GenMuon"+std::to_string(i+1)+"phi/F").c_str());
    }    
    
         
    // reco objects
    for(int i = 0; i < WriteNElectrons; i++) tree->Branch(("Electron"+std::to_string(i+1)).c_str(), &(Electrons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNMuons; i++) tree->Branch(("Muon"+std::to_string(i+1)).c_str(), &(Muons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNLeptons; i++) tree->Branch(("Lepton"+std::to_string(i+1)).c_str(), &(Leptons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNTaus; i++) tree->Branch(("Tau"+std::to_string(i+1)).c_str(), &(Taus[i]), ObjectsFormat::ListTauType().c_str());
    for(int i = 0; i < WriteNPhotons; i++) tree->Branch(("Photon"+std::to_string(i+1)).c_str(), &(Photons[i]), ObjectsFormat::ListPhotonType().c_str());
    for(int i = 0; i < WriteNJets; i++) tree->Branch(("Jet"+std::to_string(i+1)).c_str(), &(Jets[i]), ObjectsFormat::ListJetType().c_str());
    tree->Branch("MEt", &MEt, ObjectsFormat::ListMEtType().c_str());
    
}

// ------------ method called once each job just after ending the event loop  ------------
void SSleptons::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SSleptons::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SSleptons);
