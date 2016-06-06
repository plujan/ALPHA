// -*- C++ -*-
//
// Package:    Analysis/Diboson
// Class:      Diboson
// 
/**\class Diboson Diboson.cc Analysis/Diboson/plugins/Diboson.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alberto Zucchetta
//         Created:  Thu, 28 Apr 2016 08:28:54 GMT
//
//

#include "Diboson.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Diboson::Diboson(const edm::ParameterSet& iConfig):
    GenPSet(iConfig.getParameter<edm::ParameterSet>("genSet")),
    PileupPSet(iConfig.getParameter<edm::ParameterSet>("pileupSet")),
    TriggerPSet(iConfig.getParameter<edm::ParameterSet>("triggerSet")),
    ElectronPSet(iConfig.getParameter<edm::ParameterSet>("electronSet")),
    MuonPSet(iConfig.getParameter<edm::ParameterSet>("muonSet")),
    TauPSet(iConfig.getParameter<edm::ParameterSet>("tauSet")),
    PhotonPSet(iConfig.getParameter<edm::ParameterSet>("photonSet")),
    JetPSet(iConfig.getParameter<edm::ParameterSet>("jetSet")),
    FatJetPSet(iConfig.getParameter<edm::ParameterSet>("fatJetSet")),
    WriteNElectrons(iConfig.getParameter<int>("writeNElectrons")),
    WriteNMuons(iConfig.getParameter<int>("writeNMuons")),
    WriteNLeptons(iConfig.getParameter<int>("writeNLeptons")),
    WriteNTaus(iConfig.getParameter<int>("writeNTaus")),
    WriteNPhotons(iConfig.getParameter<int>("writeNPhotons")),
    WriteNJets(iConfig.getParameter<int>("writeNJets")),
    WriteNFatJets(iConfig.getParameter<int>("writeNFatJets")),
    HistFile(iConfig.getParameter<std::string>("histFile")),
    Verbose(iConfig.getParameter<bool>("verbose"))
{
    //now do what ever initialization is needed
    usesResource("TFileService");
    
    // Initialize Objects
    theGenAnalyzer      = new GenAnalyzer(GenPSet, consumesCollector());
    thePileupAnalyzer   = new PileupAnalyzer(PileupPSet, consumesCollector());
    theTriggerAnalyzer  = new TriggerAnalyzer(TriggerPSet, consumesCollector());
    theElectronAnalyzer = new ElectronAnalyzer(ElectronPSet, consumesCollector());
    theMuonAnalyzer     = new MuonAnalyzer(MuonPSet, consumesCollector());
    theTauAnalyzer      = new TauAnalyzer(TauPSet, consumesCollector());
    thePhotonAnalyzer   = new PhotonAnalyzer(PhotonPSet, consumesCollector());
    theJetAnalyzer      = new JetAnalyzer(JetPSet, consumesCollector());
    theFatJetAnalyzer   = new JetAnalyzer(FatJetPSet, consumesCollector());
    //theBTagAnalyzer     = new BTagAnalyzer(BTagAlgo);
    
    std::vector<std::string> TriggerList(TriggerPSet.getParameter<std::vector<std::string> >("paths"));
    for(unsigned int i = 0; i < TriggerList.size(); i++) TriggerMap[ TriggerList[i] ] = false;
    
    
    // ---------- Plots Initialization ----------
    TFileDirectory allDir=fs->mkdir("All/");
    TFileDirectory genDir=fs->mkdir("Gen/");
    TFileDirectory eleDir=fs->mkdir("Electrons/");
    TFileDirectory muoDir=fs->mkdir("Muons/");
    TFileDirectory jetDir=fs->mkdir("Jets/");
    TFileDirectory kinDir=fs->mkdir("Kin/");
    
    // Make TH1F
    //    std::vector<std::string> nLabels={"Trigger", "Lep #geq 2", "Z cand ", "Jets #geq 2", "Z mass ", "bJets #geq 1", "bJets #geq 2", "h mass ", "#slash{E}_{T}", "Final"};
    std::vector<std::string> nLabels={"All", "Z Lep", "H Merged", "H Resolved", "X Merged", "X Resolved", "..", "..", "..", ".."};
    
    int nbins;
    float min, max;
    std::string name, title, opt;
    
    ifstream histFile(HistFile);
    if(!histFile.is_open()) {
        throw cms::Exception("Diboson Analyzer", HistFile + " file not found");
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


Diboson::~Diboson() {
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
    delete theFatJetAnalyzer;
    //delete theBTagAnalyzer;
    
}


//
// member functions
//

// ------------ method called for each event  ------------
void Diboson::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    isMC = !iEvent.isRealData();
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
    for(int i = 0; i < WriteNFatJets; i++) ObjectsFormat::ResetFatJetType(FatJets[i]);
    ObjectsFormat::ResetMEtType(MEt);
    ObjectsFormat::ResetCandidateType(V);
    ObjectsFormat::ResetCandidateType(HMerged);
    ObjectsFormat::ResetCandidateType(HResolved);
    
    
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
    theJetAnalyzer->CleanJetsFromMuons(JetsVect, MuonVect, 0.4);
    theJetAnalyzer->CleanJetsFromElectrons(JetsVect, ElecVect, 0.4);
    // Fat Jets
    std::vector<pat::Jet> FatJetsVect = theFatJetAnalyzer->FillJetVector(iEvent);
    theFatJetAnalyzer->CleanJetsFromMuons(FatJetsVect, MuonVect, 1.);
    theFatJetAnalyzer->CleanJetsFromElectrons(FatJetsVect, ElecVect, 1.);
    // Missing Energy
    pat::MET MET = theJetAnalyzer->FillMetVector(iEvent);
    
    // Gen weights
    std::map<std::string, float> GenWeight = theGenAnalyzer->FillWeightsMap(iEvent);
    EventWeight *= GenWeight["event"];
    // Gen Particles
    std::vector<reco::GenParticle> GenPVect = theGenAnalyzer->FillGenVector(iEvent);
    // Gen candidates
    //reco::Candidate* theGenZ = theGenAnalyzer->FindGenParticle(GenPVect, 23);
    std::vector<int> LepIds = {11,13};
    std::vector<int> HadIds = {1,2,3,4,5};
    reco::GenParticle* theGenLep = theGenAnalyzer->FindGenParticleGenByIds(GenPVect, LepIds);
    reco::GenParticle* theGenHad = theGenAnalyzer->FindGenParticleGenByIds(GenPVect, HadIds);

    //Gen level plots and candidates
    double GenZLepMass = 0.;
    double GenZHadMass = 0.;
    double GenXMass = 0.;
    bool isGenZZ = false;
    if(theGenLep!=NULL && theGenHad!=NULL){
        const reco::Candidate* theGenZLep = theGenAnalyzer->FindMother(theGenLep);
        const reco::Candidate* theGenZHad = theGenAnalyzer->FindMother(theGenHad);
        if(theGenZLep!=NULL && theGenZLep->pdgId()==23 && theGenZHad!=NULL && theGenZHad->pdgId()==23){
            Hist["g_ZLepMass"]->Fill(theGenZLep->mass(), EventWeight);
            Hist["g_ZLepPt"]->Fill(theGenZLep->pt(), EventWeight);
            Hist["g_LepPt"]->Fill(theGenLep->pt(), EventWeight);
            GenZLepMass = theGenZLep->mass();
            Hist["g_ZHadMass"]->Fill(theGenZHad->mass(), EventWeight);
            Hist["g_ZHadPt"]->Fill(theGenZHad->pt(), EventWeight);
            Hist["g_HadPt"]->Fill(theGenHad->pt(), EventWeight);
            Hist["g_HadEta"]->Fill(theGenHad->eta(), EventWeight);
            Hist["g_ZZDR"]->Fill(reco::deltaR(theGenZHad->eta(),theGenZHad->phi(),theGenZLep->eta(),theGenZLep->phi()), EventWeight);
            Hist["g_ZZDPhi"]->Fill(reco::deltaPhi(theGenZHad->phi(),theGenZLep->phi()), EventWeight);
            Hist["g_LepHadDR"]->Fill(reco::deltaR(theGenHad->eta(),theGenHad->phi(),theGenLep->eta(),theGenLep->phi()), EventWeight);
            GenZHadMass = theGenZHad->mass();
            isGenZZ = true;
        }
    }

    
    reco::Candidate* theGenX = theGenAnalyzer->FindGenParticleByIdAndStatus(GenPVect, 39, 62);
    if(theGenX!=NULL && theGenLep!=NULL && theGenHad!=NULL && isGenZZ){
        Hist["g_XMass"]->Fill(theGenX->mass(), EventWeight);
        Hist["g_XPt"]->Fill(theGenX->pt(), EventWeight);
        GenXMass = theGenX->mass();
    }
    

    // Pu weight
    PUWeight = thePileupAnalyzer->GetPUWeight(iEvent);
    EventWeight *= PUWeight;
    
    // Trigger
    theTriggerAnalyzer->FillTriggerMap(iEvent, TriggerMap);
    EventWeight *= TriggerWeight;
    
    
    Hist["a_nEvents"]->Fill(1., EventWeight);
    Hist["e_nEvents"]->Fill(1., EventWeight);
    Hist["m_nEvents"]->Fill(1., EventWeight);
    
    // ---------- Do analysis selections ----------
    
    // ---------- Z TO LEPTONS ----------
    isZtoEE = isZtoMM = isWtoEN = isWtoMN = isZtoNN = false;
    
    if(MuonVect.size()>=2 && ElecVect.size()>=2) {
        if(MuonVect.at(0).pt() > ElecVect.at(0).pt()) isZtoMM=true;
        else isZtoEE=true;
    }
    else if(ElecVect.size()>=2) isZtoEE=true;
    else if(MuonVect.size()>=2) isZtoMM=true;
    //else {if(Verbose) std::cout << " - No Iso SF OS Leptons" << std::endl; return;}

    // ---------- W TO LEPTON and NEUTRINO ----------
    
    else if(MuonVect.size()==1 && ElecVect.size()==1) {
        if(MuonVect.at(0).pt() > ElecVect.at(0).pt()) isWtoMN=true;
        else isWtoEN=true;
    }
    else if(ElecVect.size()==1) isWtoEN=true;
    else if(MuonVect.size()==1) isWtoMN=true;
    //else {if(Verbose) std::cout << " - No W to Leptons" << std::endl; return;}

    // ----------- Z TO NEUTRINOS ---------------

    else { if(Verbose) std::cout << " - No charged leptons" << std::endl;}// return;}

    // AZh lepton choice
//    int l1(0), l2(-1);
//    if(isZtoEE) {
//        for(unsigned int i=1; i<ElecVect.size(); i++) if(l2<0 && ElecVect.at(i).charge()!=ElecVect.at(l1).charge()) l2=i;
//    }
//    else if(isZtoMM) {
//        for(unsigned int i=1; i<MuonVect.size(); i++) if(l2<0 && MuonVect.at(i).charge()!=MuonVect.at(l1).charge()) l2=i;
//    }
//    if(l1<0 || l2<0) {if(Verbose) std::cout << " - No OS SF leptons" << std::endl; return;}

    // Reconstruct V candidate
    pat::CompositeCandidate theV;
    if(isZtoMM) {
        if(MuonVect.at(0).charge()*MuonVect.at(1).charge()<0){
            theV.addDaughter(MuonVect.at(0));
            theV.addDaughter(MuonVect.at(1));
            addP4.set(theV);
            Hist["a_nEvents"]->Fill(2., EventWeight);
            Hist["m_nEvents"]->Fill(2., EventWeight);
//            MuonVect.at(0).addUserFloat("trkIso",theMuonAnalyzer->FixTrackerIsolation(MuonVect.at(0),MuonVect.at(1)).at(0));
//            MuonVect.at(1).addUserFloat("trkIso",theMuonAnalyzer->FixTrackerIsolation(MuonVect.at(0),MuonVect.at(1)).at(1));
        }
        else { if(Verbose) std::cout << " - No OS SF leptons" << std::endl; return;  }
    }
    else if(isZtoEE){
        if(ElecVect.at(0).charge()*ElecVect.at(1).charge()<0){
            theV.addDaughter(ElecVect.at(0));
            theV.addDaughter(ElecVect.at(1));
            addP4.set(theV);
            Hist["a_nEvents"]->Fill(2., EventWeight);
            Hist["e_nEvents"]->Fill(2., EventWeight);
        }
        else { if(Verbose) std::cout << " - No OS SF leptons" << std::endl; return;  }
    }
    else{
        if(Verbose) std::cout << "No dilepton!" << std::endl;
        return;
    }
    /*
    else if(isWtoMN){
        theV = createkW(MuonVect.at(0), MET);
    }
    else if(isWtoEN){
        theV = createkW(ElecVect.at(0), MET);
    }
    else {
        //Z to nu nu
    }
    */
    // ---------- Z TO HADRONS ----------
    pat::CompositeCandidate theHMerged;
    pat::CompositeCandidate theHResolved;
    pat::CompositeCandidate theH;
    isMerged = isResolved = false;
    /////////////////// Highest pT method ////////////////////
    
    // Resolved topology
    if(JetsVect.size() < 2) {if(Verbose) std::cout << " - N jets < 2" << std::endl;} // return;}
    else {
        theH.addDaughter(JetsVect.at(0));
        theH.addDaughter(JetsVect.at(1));
        addP4.set(theH);
  //std::cout << "resolved mass: " << theH.mass() << std::endl;
        //if(theH.mass()<40 || theH.pt()<100) theH.clearDaughters();
        //Hist["a_HAK4Mass_HPt"]->Fill(theH.mass(), EventWeight);
    }

    // Boosted topology
    if(FatJetsVect.size() < 1) {if(Verbose) std::cout << " - N fat jets < 1" << std::endl;}
    else {
      //Hist["a_HAK8Mass_HPt"]->Fill(FatJetsVect.at(0).hasUserFloat("ak8PFJetsCHSSoftDropMass")?FatJetsVect.at(0).userFloat("ak8PFJetsCHSSoftDropMass"):FatJetsVect.at(0).mass(), EventWeight);
      //std::cout << "merged mass: " << FatJetsVect.at(0).mass() << std::endl;
      if(theH.pt()<FatJetsVect.at(0).pt()){
        theH.clearDaughters();
        theH.addDaughter(FatJetsVect.at(0));
        addP4.set(theH);
        theH.addUserFloat("softdropMass",FatJetsVect.at(0).userFloat("ak8PFJetsCHSSoftDropMass"));
        //if(theH.mass()>180) theH.clearDaughters();
      }
    }
    Hist["a_HMass_HPt"]->Fill(theH.hasUserFloat("softdropMass") ? theH.userFloat("softdropMass") : theH.mass(), EventWeight);
    //std::cout << "chosen mass: " << theH.mass() << std::endl;
    
    // Reset theH
    //theH.clearDaughters();
    

    /////////////////// Prefer merged AK8 jet method ////////////////////
    
    // Boosted topology
    //if(FatJetsVect.size() < 1) {if(Verbose) std::cout << " - N fat jets < 1" << std::endl;}

    if(FatJetsVect.size()>=1){
      isMerged = true;
      theHMerged.addDaughter(FatJetsVect.at(0));
      theHMerged.addUserFloat("softdropMass",FatJetsVect.at(0).userFloat("ak8PFJetsCHSSoftDropMass"));
      addP4.set(theHMerged);
      Hist["a_nEvents"]->Fill(3., EventWeight);
      if(isZtoEE) Hist["e_nEvents"]->Fill(3., EventWeight);
      if(isZtoMM)Hist["m_nEvents"]->Fill(3., EventWeight);
      //std::cout << "merged mass: " << theH.mass() << std::endl;
      Hist["a_HAK8Mass"]->Fill(FatJetsVect.at(0).hasUserFloat("ak8PFJetsCHSSoftDropMass")?FatJetsVect.at(0).userFloat("ak8PFJetsCHSSoftDropMass"):FatJetsVect.at(0).mass(), EventWeight);
    }

    //else{
        // Resolved topology if we don't have the right AK8
        //if(JetsVect.size() < 2) {if(Verbose) std::cout << " - N fat jet <1 or with pruned mass < 30 GeV and N ak4 jets < 2, no H candidate" << std::endl;} // return;}
    if(JetsVect.size() >= 2){
            isResolved = true;
            theHResolved.addDaughter(JetsVect.at(0));
            theHResolved.addDaughter(JetsVect.at(1));
            addP4.set(theHResolved);
            Hist["a_nEvents"]->Fill(4., EventWeight);
            if(isZtoEE)Hist["e_nEvents"]->Fill(4., EventWeight);
            if(isZtoMM)Hist["m_nEvents"]->Fill(4., EventWeight);
      //std::cout << "resolved mass: " << theH.mass() << std::endl;
            Hist["a_HAK4Mass"]->Fill(theHResolved.mass(), EventWeight);
        }
      //}

    if(FatJetsVect.size()>=1) Hist["a_HMass_PM"]->Fill(theHMerged.userFloat("softdropMass"), EventWeight);
    else if(FatJetsVect.size()<1 && JetsVect.size() >= 2) Hist["a_HMass_PM"]->Fill(theHResolved.mass(), EventWeight);
    //std::cout << "chosen mass: " << theH.mass() << std::endl;
    

    // Global candidate
    pat::CompositeCandidate theXMerged;
    pat::CompositeCandidate theXResolved;
    if((isZtoEE || isZtoMM) && theHMerged.numberOfDaughters()>0){
        theXMerged.addDaughter(theHMerged);
        theXMerged.addDaughter(theV);
        addP4.set(theXMerged);
        Hist["a_nEvents"]->Fill(5., EventWeight);
        if(isZtoMM) Hist["m_nEvents"]->Fill(5., EventWeight);
        if(isZtoEE) Hist["e_nEvents"]->Fill(5., EventWeight);
   }

    if((isZtoEE || isZtoMM) && theHResolved.numberOfDaughters()>0){
        theXResolved.addDaughter(theHResolved);
        theXResolved.addDaughter(theV);
        addP4.set(theXResolved);
        Hist["a_nEvents"]->Fill(6., EventWeight);
        if(isZtoMM) Hist["m_nEvents"]->Fill(6., EventWeight);
        if(isZtoEE) Hist["e_nEvents"]->Fill(6., EventWeight);
    }

    else{
        if(Verbose) std::cout << "No dilepton, no X!" << std::endl;
        return;
    }

    /*
    else if(theH.numberOfDaughters()>0){//if is Z to nu nu: apply recoil mass formula
        theX = recoilMassFormula(theH,MET);
    }

    else{//no H, no X
        if(Verbose) std::cout << "No H, no X" << std::endl;
        return;
    }
    */

    // ---------- Print Summary ----------
    if(Verbose) {
        std::cout << " --- Event n. " << iEvent.id().event() << ", lumi " << iEvent.luminosityBlock() << ", run " << iEvent.id().run() << ", weight " << EventWeight << std::endl;
        std::cout << "number of electrons: " << ElecVect.size() << std::endl;
        for(unsigned int i = 0; i < ElecVect.size(); i++) std::cout << "  electron [" << i << "]\tpt: " << ElecVect[i].pt() << "\teta: " << ElecVect[i].eta() << "\tphi: " << ElecVect[i].phi() << "\tmass: " << ElecVect[i].mass() << "\tcharge: " << ElecVect[i].charge() << std::endl;
        std::cout << "number of muons:     " << MuonVect.size() << std::endl;
        for(unsigned int i = 0; i < MuonVect.size(); i++) std::cout << "  muon     [" << i << "]\tpt: " << MuonVect[i].pt() << "\teta: " << MuonVect[i].eta() << "\tphi: " << MuonVect[i].phi() << "\tmass: " << MuonVect[i].mass() << "\tcharge: " << MuonVect[i].charge() << std::endl;
        std::cout << "number of taus:  " << TauVect.size() << std::endl;
        for(unsigned int i = 0; i < TauVect.size(); i++) std::cout << "  tau  [" << i << "]\tpt: " << TauVect[i].pt() << "\teta: " << TauVect[i].eta() << "\tphi: " << TauVect[i].phi() << std::endl;
        std::cout << "number of photons:  " << PhotonVect.size() << std::endl;
        for(unsigned int i = 0; i < PhotonVect.size(); i++) std::cout << "  photon  [" << i << "]\tpt: " << PhotonVect[i].pt() << "\teta: " << PhotonVect[i].eta() << "\tphi: " << PhotonVect[i].phi() << std::endl;
        std::cout << "number of AK4 jets:  " << JetsVect.size() << std::endl;
        for(unsigned int i = 0; i < JetsVect.size(); i++) std::cout << "  AK4 jet  [" << i << "]\tpt: " << JetsVect[i].pt() << "\teta: " << JetsVect[i].eta() << "\tphi: " << JetsVect[i].phi() << "\tmass: " << JetsVect[i].mass() << std::endl;
        std::cout << "number of AK8 jets:  " << FatJetsVect.size() << std::endl;
        for(unsigned int i = 0; i < FatJetsVect.size(); i++) std::cout << "  AK8 jet  [" << i << "]\tpt: " << FatJetsVect[i].pt() << "\teta: " << FatJetsVect[i].eta() << "\tphi: " << FatJetsVect[i].phi() << "\tmass: " << FatJetsVect[i].mass() << std::endl;
        std::cout << "Missing energy:      " << MET.pt() << std::endl;
        std::cout << "V leptonic mass:     " << theV.mass() << ", generated: " << GenZLepMass << std::endl;
        std::cout << "Z merged hadronic mass:     " << theHMerged.mass() << ", generated: " << GenZHadMass << std::endl;
        std::cout << "X merged candidate mass:    " << theXMerged.mass() << ", generated: " << GenXMass << std::endl;
        std::cout << "Z resolved hadronic mass:     " << theHResolved.mass() << ", generated: " << GenZHadMass << std::endl;
        std::cout << "X resolved candidate mass:    " << theXResolved.mass() << ", generated: " << GenXMass << std::endl;
    }

    
    // ---------- Fill objects ----------
    if(isZtoEE || isWtoEN) for(unsigned int i = 0; i < Leptons.size() && i < ElecVect.size(); i++) ObjectsFormat::FillElectronType(Leptons[i], &ElecVect[i], isMC);
    else if(isZtoMM || isWtoMN) for(unsigned int i = 0; i < Leptons.size() && i < MuonVect.size(); i++) ObjectsFormat::FillMuonType(Leptons[i], &MuonVect[i], isMC);
    for(unsigned int i = 0; i < Taus.size() && i < TauVect.size(); i++) ObjectsFormat::FillTauType(Taus[i], &TauVect[i], isMC);
    for(unsigned int i = 0; i < Photons.size() && i < PhotonVect.size(); i++) ObjectsFormat::FillPhotonType(Photons[i], &PhotonVect[i], isMC);
    for(unsigned int i = 0; i < Jets.size() && i < JetsVect.size(); i++) ObjectsFormat::FillJetType(Jets[i], &JetsVect[i], isMC);
    for(unsigned int i = 0; i < FatJets.size() && i < FatJetsVect.size(); i++) ObjectsFormat::FillFatJetType(FatJets[i], &FatJetsVect[i], isMC);
    ObjectsFormat::FillMEtType(MEt, &MET, isMC);
    ObjectsFormat::FillCandidateType(V, &theV, isMC);
    ObjectsFormat::FillCandidateType(HMerged, &theHMerged, isMC);
    ObjectsFormat::FillCandidateType(XMerged, &theXMerged, isMC);
    ObjectsFormat::FillCandidateType(HResolved, &theHResolved, isMC);
    ObjectsFormat::FillCandidateType(XResolved, &theXResolved, isMC);
        
    // Lepton and Trigger SF
    if(isMC) {
          if(isZtoEE) {
            //TriggerWeight*=theElectronAnalyzer->GetDoubleElectronTriggerSF(ElecVect.at(0), ElecVect.at(1));
            LeptonWeight*=theElectronAnalyzer->GetElectronIdSF(ElecVect.at(0), ElectronPSet.getParameter<int>("electron1id"));
            LeptonWeight*=theElectronAnalyzer->GetElectronIdSF(ElecVect.at(1), ElectronPSet.getParameter<int>("electron2id"));
            LeptonWeight*=theElectronAnalyzer->GetElectronRecoEffSF(ElecVect.at(0));
            LeptonWeight*=theElectronAnalyzer->GetElectronRecoEffSF(ElecVect.at(1));
      }
        if(isZtoMM) {
      //TriggerWeight*=theMuonAnalyzer->GetDoubleMuonTriggerSF(MuonVect.at(0), MuonVect.at(1));
            LeptonWeight*=theMuonAnalyzer->GetMuonIdSF(MuonVect.at(0), MuonPSet.getParameter<int>("muon1id"));
            LeptonWeight*=theMuonAnalyzer->GetMuonIdSF(MuonVect.at(1), MuonPSet.getParameter<int>("muon2id"));
            LeptonWeight*=theMuonAnalyzer->GetMuonIsoSF(MuonVect.at(0), MuonPSet.getParameter<int>("muon1iso"));
            LeptonWeight*=theMuonAnalyzer->GetMuonIsoSF(MuonVect.at(1), MuonPSet.getParameter<int>("muon2iso"));
        }
        if(isWtoEN) {
            //TriggerWeight*=theElectronAnalyzer->GetDoubleElectronTriggerSF(ElecVect.at(0));
      LeptonWeight*=theElectronAnalyzer->GetElectronIdSF(ElecVect.at(0), ElectronPSet.getParameter<int>("electron1id"));
            LeptonWeight*=theElectronAnalyzer->GetElectronRecoEffSF(ElecVect.at(0));
      }
        if(isWtoMN) {
      //TriggerWeight*=theMuonAnalyzer->GetDoubleMuonTriggerSF(MuonVect.at(0));
            LeptonWeight*=theMuonAnalyzer->GetMuonIdSF(MuonVect.at(0), MuonPSet.getParameter<int>("muon1iso"));
            LeptonWeight*=theMuonAnalyzer->GetMuonIsoSF(MuonVect.at(0), MuonPSet.getParameter<int>("muon1iso"));
        }
    }
    //EventWeight*=TriggerWeight;
    EventWeight*=LeptonWeight;

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
//        std::cout << "\tReconstructed Z candidate from " << (isZtoMM ? "muons" : "electrons") << " " << l1 << " and " << l2 << " with mass: " << theV.mass() << std::endl;
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
//    theJetAnalyzer->ApplyRecoilCorrections(MET, &MET.genMET()->p4(), &theV.p4(), 0);

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
void Diboson::beginJob() {
    
    // Object objects are created only one in the begin job. The reference passed to the branch has to be the same
    for(int i = 0; i < WriteNElectrons; i++) Electrons.push_back( LeptonType() );
    for(int i = 0; i < WriteNMuons; i++) Muons.push_back( LeptonType() );
    for(int i = 0; i < WriteNLeptons; i++) Leptons.push_back( LeptonType() );
    for(int i = 0; i < WriteNTaus; i++) Taus.push_back( TauType() );
    for(int i = 0; i < WriteNPhotons; i++) Photons.push_back( PhotonType() );
    for(int i = 0; i < WriteNJets; i++) Jets.push_back( JetType() );
    for(int i = 0; i < WriteNFatJets; i++) FatJets.push_back( FatJetType() );
    
    
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
    
    // Set trigger branches
    for(auto it = TriggerMap.begin(); it != TriggerMap.end(); it++) tree->Branch(it->first.c_str(), &(it->second), (it->first+"/O").c_str());
    
    // Analysis variables
    tree->Branch("isZtoEE", &isZtoEE, "isZtoEE/O");
    tree->Branch("isZtoMM", &isZtoMM, "isZtoMM/O");
    tree->Branch("isWtoEN", &isWtoEN, "isWtoEN/O");
    tree->Branch("isWtoMN", &isWtoMN, "isWtoMN/O");
    tree->Branch("isZtoNN", &isZtoNN, "isZtoNN/O");
    tree->Branch("isMerged", &isMerged, "isMerged/O");
    tree->Branch("isResolved", &isResolved, "isResolved/O");
    
    // Set Branches for objects
    for(int i = 0; i < WriteNElectrons; i++) tree->Branch(("Electron"+std::to_string(i+1)).c_str(), &(Electrons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNMuons; i++) tree->Branch(("Muon"+std::to_string(i+1)).c_str(), &(Muons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNLeptons; i++) tree->Branch(("Lepton"+std::to_string(i+1)).c_str(), &(Leptons[i]), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNTaus; i++) tree->Branch(("Tau"+std::to_string(i+1)).c_str(), &(Taus[i]), ObjectsFormat::ListTauType().c_str());
    for(int i = 0; i < WriteNPhotons; i++) tree->Branch(("Photon"+std::to_string(i+1)).c_str(), &(Photons[i]), ObjectsFormat::ListPhotonType().c_str());
    for(int i = 0; i < WriteNJets; i++) tree->Branch(("Jet"+std::to_string(i+1)).c_str(), &(Jets[i]), ObjectsFormat::ListJetType().c_str());
    for(int i = 0; i < WriteNFatJets; i++) tree->Branch(("FatJet"+std::to_string(i+1)).c_str(), &(FatJets[i]), ObjectsFormat::ListFatJetType().c_str());
    tree->Branch("MEt", &MEt, ObjectsFormat::ListMEtType().c_str());
    tree->Branch("V", &V, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("HMerged", &HMerged, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("XMerged", &XMerged, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("HResolved", &HResolved, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("XResolved", &XResolved, ObjectsFormat::ListCandidateType().c_str());
}

// ------------ method called once each job just after ending the event loop  ------------
void Diboson::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Diboson::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

pat::CompositeCandidate Diboson::createkW(reco::Candidate& lep, pat::MET& met){
    pat::CompositeCandidate thekW;
    // W kinematical reconstruction
    float pz = 0.;
    float a = pow(80.4,2) - pow(lep.mass(),2) + 2.*lep.px()*met.px() + 2.*lep.py()*met.py();
    float A = 4*( pow(lep.energy(),2) - pow(lep.pz(),2) );
    float B = -4*a*lep.pz();
    float C = 4*pow(lep.energy(),2) * (pow(met.px(),2)  + pow(met.py(),2)) - pow(a,2);
    float D = pow(B,2) - 4*A*C;
    // If there are real solutions, use the one with lowest pz                                            
    if (D>=0){
        float s1 = (-B+sqrt(D))/(2*A);
        float s2 = (-B-sqrt(D))/(2*A);
        if(fabs(s1)<fabs(s2)) pz=s1;
        else pz=s2;
    }
    // Otherwise, use real part                                                                           
    else{
         pz = -B/(2*A);
    }
    reco::Particle::LorentzVector neutrino( met.px(), met.py(), pz, sqrt( pow(met.pt(),2)+pow(pz,2) ) );
    reco::Particle::LorentzVector kW = lep.p4() + neutrino;
    thekW.setP4(kW);
    thekW.setCharge(lep.charge());
    return thekW;
}

pat::CompositeCandidate Diboson::recoilMassFormula(pat::CompositeCandidate& H, pat::MET& met){
    pat::CompositeCandidate X;
    AddFourMomenta addP4;
    X.addDaughter(H);
    X.addDaughter(met);
    addP4.set(X);
    reco::Particle::LorentzVector metp4 = met.p4();
    reco::Particle::LorentzVector Xp4;
    metp4.SetPz(-H.pz());
    Xp4 += metp4;
    Xp4.SetPz(0);
    float B = -2.*H.energy();
    float C = pow(H.mass(),2) - pow(90.18,2);
    float D = pow(B,2) - 4*1*C;
    float mX;
    if(D>0){
        float s1 = (-B+sqrt(D))/2.;
        float s2 = (-B-sqrt(D))/2.;
        if(fabs(s1)>fabs(s2)) mX = s1;
        else mX = s2;
    }
    else{
        mX = -B/2.;
    }
    Xp4.SetE(sqrt(pow(mX,2) + pow(Xp4.Px(),2) + pow(Xp4.Py(),2) + pow(Xp4.Pz(),2)));
    X.setP4(Xp4);
    X.setCharge(0);
    return X;
}


//define this as a plug-in
DEFINE_FWK_MODULE(Diboson);
