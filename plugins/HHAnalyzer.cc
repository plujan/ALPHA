// -*- C++ -*-
//
// Package:    Analysis/HHAnalyzer
// Class:      HHAnalyzer
// 
/**\class HHAnalyzer HHAnalyzer.cc Analysis/HHAnalyzer/plugins/HHAnalyzer.cc
 Description: [one line class summary]
 Implementation:     [Notes on implementation]
*/
//
// Original Author:  Alberto Zucchetta
// Created:  Thu, 28 Apr 2016 08:28:54 GMT
// Modified for HH analysis - Sept 2016
//

#include "HHAnalyzer.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HHAnalyzer::HHAnalyzer(const edm::ParameterSet& iConfig):
    GenPSet(iConfig.getParameter<edm::ParameterSet>("genSet")),
    PileupPSet(iConfig.getParameter<edm::ParameterSet>("pileupSet")),
    TriggerPSet(iConfig.getParameter<edm::ParameterSet>("triggerSet")),
    ElectronPSet(iConfig.getParameter<edm::ParameterSet>("electronSet")),
    MuonPSet(iConfig.getParameter<edm::ParameterSet>("muonSet")),
    JetPSet(iConfig.getParameter<edm::ParameterSet>("jetSet")),
    FatJetPSet(iConfig.getParameter<edm::ParameterSet>("fatJetSet")),
    WriteNElectrons(iConfig.getParameter<int>("writeNElectrons")),
    WriteNMuons(iConfig.getParameter<int>("writeNMuons")),
    WriteNLeptons(iConfig.getParameter<int>("writeNLeptons")),
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
    std::vector<std::string> nLabels={"All", "Trigger", "Iso Lep #geq 2", "Z cand ", "Jets #geq 2", "Z mass ", "h mass ", "Top veto", "bJets #geq 1", "bJets #geq 2"};
    
    int nbins;
    float min, max;
    std::string name, title, opt;
    
    ifstream histFile(HistFile);
    if(!histFile.is_open()) {
        throw cms::Exception("HHAnalyzer Analyzer", HistFile + " file not found");
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


HHAnalyzer::~HHAnalyzer() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    std::cout << "---------- ENDING  ----------" << std::endl;
    
    delete theGenAnalyzer;
    delete thePileupAnalyzer;
    delete theTriggerAnalyzer;
    delete theElectronAnalyzer;
    delete theMuonAnalyzer;
    delete theJetAnalyzer;
    delete theFatJetAnalyzer;
    //delete theBTagAnalyzer;    
}


//
// member functions
//

// ------------ method called for each event  ------------
void HHAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    isMC = !iEvent.isRealData();
    EventNumber = iEvent.id().event();
    LumiNumber = iEvent.luminosityBlock();
    RunNumber = iEvent.id().run();
    
    EventWeight = StitchWeight = ZewkWeight = WewkWeight = PUWeight = TriggerWeight = LeptonWeight = 1.;
    FacWeightUp = FacWeightDown = RenWeightUp = RenWeightDown = ScaleWeightUp = ScaleWeightDown = 1.;
    isZtoEE = isZtoMM = isTtoEM = isWtoEN = isWtoMN = isZtoNN = false;
    nPV = nElectrons = nMuons = nJets = nFatJets = nBTagJets = 1;
    nVetoElectrons = nLooseMuons = 0;
    MaxJetBTag = MaxFatJetBTag = Chi2 = -1.;
    MinJetMetDPhi = 10.;
    massRecoilFormula = -1.;
    /*
    Lepton1_isMuon = Lepton1_isElectron = Lepton1_isLoose = Lepton1_isHighPt = Lepton1_isTrackerHighPt = Lepton1_isTight = false;
    Lepton1_pt = Lepton1_trkIso = -1.;
    Lepton2_isMuon = Lepton2_isElectron = Lepton2_isLoose = Lepton2_isHighPt = Lepton2_isTrackerHighPt = Lepton2_isTight = false;
    Lepton2_pt = Lepton2_trkIso = -1.;
    MEt_pt = -1.;
    V_pt = V_dPhi = V_mass = V_tmass = -1.;
    X_pt = X_dPhi = X_mass = X_tmass = -1.;
    FatJet1_isTight = false;
    FatJet1_pt = FatJet1_prunedMass = FatJet1_softdropMass = FatJet1_softdropPuppiMass = FatJet1_prunedMassCorr = FatJet1_softdropMassCorr = FatJet1_softdropPuppiMassCorr = FatJet1_chsTau21 = FatJet1_puppiTau21 = FatJet1_ddtTau21 = FatJet1_CSV1 = FatJet1_CSV2 = -1.;
    */
    AddFourMomenta addP4;
    
    // Initialize types
    for(int i = 0; i < WriteNElectrons; i++) ObjectsFormat::ResetLeptonType(Electrons[i]);
    for(int i = 0; i < WriteNMuons; i++) ObjectsFormat::ResetLeptonType(Muons[i]);
    for(int i = 0; i < WriteNLeptons; i++) ObjectsFormat::ResetLeptonType(Leptons[i]);
    for(int i = 0; i < WriteNFatJets; i++) ObjectsFormat::ResetFatJetType(FatJets[i]);
    ObjectsFormat::ResetMEtType(MEt);
    
    ObjectsFormat::ResetCandidateType(V);
    ObjectsFormat::ResetCandidateType(X);
    /*ObjectsFormat::ResetCandidateType(H);
    ObjectsFormat::ResetCandidateType(A);
    ObjectsFormat::ResetCandidateType(HMerged);
    ObjectsFormat::ResetCandidateType(XMerged);
    ObjectsFormat::ResetCandidateType(HResolved);
    ObjectsFormat::ResetCandidateType(XResolved);
    ObjectsFormat::ResetCandidateType(HResolvedPt);
    ObjectsFormat::ResetCandidateType(XResolvedPt);
    ObjectsFormat::ResetCandidateType(HResolvedHpt);
    ObjectsFormat::ResetCandidateType(XResolvedHpt);
    ObjectsFormat::ResetCandidateType(HResolvedDZ);
    ObjectsFormat::ResetCandidateType(XResolvedDZ);
    ObjectsFormat::ResetCandidateType(HResolvedDR);
    ObjectsFormat::ResetCandidateType(XResolvedDR);
    ObjectsFormat::ResetLorentzType(kH);
    ObjectsFormat::ResetLorentzType(kA);*/
    
    Hist["a_nEvents"]->Fill(1., EventWeight);
    Hist["e_nEvents"]->Fill(1., EventWeight);
    Hist["m_nEvents"]->Fill(1., EventWeight);
    
    // -----------------------------------
    //           READ OBJECTS
    // -----------------------------------
    
    // Pu weight
    PUWeight = thePileupAnalyzer->GetPUWeight(iEvent);
    nPV = thePileupAnalyzer->GetPV(iEvent);
    Hist["a_nPVNoWeight"]->Fill(nPV, EventWeight);
    EventWeight *= PUWeight;
    Hist["a_nPVReWeight"]->Fill(nPV, EventWeight);
    
    // Trigger
    theTriggerAnalyzer->FillTriggerMap(iEvent, TriggerMap);
    EventWeight *= TriggerWeight;
    
    // Electrons
    std::vector<pat::Electron> ElecVect = theElectronAnalyzer->FillElectronVector(iEvent);
    nElectrons = ElecVect.size();
    for(unsigned int i =0; i<ElecVect.size(); i++){
        if(ElecVect.at(i).userInt("isVeto")==1) nVetoElectrons++;
    }
    // Muons
    std::vector<pat::Muon> MuonVect = theMuonAnalyzer->FillMuonVector(iEvent);
    nMuons = MuonVect.size();
    std::vector<pat::Muon> LooseMuonVect;
    for(unsigned int i =0; i<MuonVect.size(); i++){
        if(MuonVect.at(i).userInt("isLoose")==1){
	    LooseMuonVect.push_back(MuonVect.at(i));
            nLooseMuons++;
	}
    }
    // Jets
    std::vector<pat::Jet> JetsVect = theJetAnalyzer->FillJetVector(iEvent);
    theJetAnalyzer->CleanJetsFromMuons(JetsVect, MuonVect, 0.4);
    theJetAnalyzer->CleanJetsFromElectrons(JetsVect, ElecVect, 0.4);
    nJets = JetsVect.size();
    nBTagJets = theJetAnalyzer->GetNBJets(JetsVect);
    // Fat Jets
    std::vector<pat::Jet> FatJetsVect = theFatJetAnalyzer->FillJetVector(iEvent);
    //theFatJetAnalyzer->CleanJetsFromMuons(FatJetsVect, MuonVect, 1.); // Do NOT clean the fatjet now
    //theFatJetAnalyzer->CleanJetsFromElectrons(FatJetsVect, ElecVect, 1.); // Do NOT clean the fatjet now
    nFatJets = FatJetsVect.size();
    // Missing Energy
    pat::MET MET = theJetAnalyzer->FillMetVector(iEvent);
    pat::MET Neutrino(MET);
    float metNoMupx = MET.px();
    float metNoMupy = MET.py();
    for(unsigned int i=0; i<LooseMuonVect.size();i++){
      metNoMupx -= LooseMuonVect.at(i).px();
      metNoMupy -= LooseMuonVect.at(i).py();
    }
    reco::Particle::LorentzVector metNoMup4(metNoMupx, metNoMupy, 0, 0 );
    MET.addUserFloat("metNoMu",metNoMup4.px());
    MET.addUserFloat("phiNoMu",metNoMup4.phi());

    //theJetAnalyzer->ApplyRecoilCorrections(MET, &MET.genMET()->p4(), &theV.p4(), 0);
    
    
    
    // -----------------------------------
    //           GEN LEVEL
    // -----------------------------------
    
    // Gen weights
    std::map<std::string, float> GenWeight = theGenAnalyzer->FillWeightsMap(iEvent);
    EventWeight *= GenWeight["event"];
    if(GenWeight.find("2") != GenWeight.end()) FacWeightUp     = GenWeight["2"];
    if(GenWeight.find("3") != GenWeight.end()) FacWeightDown   = GenWeight["3"];
    if(GenWeight.find("4") != GenWeight.end()) RenWeightUp     = GenWeight["4"];
    if(GenWeight.find("7") != GenWeight.end()) RenWeightDown   = GenWeight["7"];
    if(GenWeight.find("5") != GenWeight.end()) ScaleWeightUp   = GenWeight["5"];
    if(GenWeight.find("9") != GenWeight.end()) ScaleWeightDown = GenWeight["9"];
    // Lhe Particles
    std::map<std::string, float> LheMap = theGenAnalyzer->FillLheMap(iEvent);
    // Mc Stitching
    StitchWeight = theGenAnalyzer->GetStitchWeight(LheMap);
    //EventWeight *= StitchWeight; // Not yet
    // Gen Particles
    std::vector<reco::GenParticle> GenPVect = theGenAnalyzer->FillGenVector(iEvent);
    // Gen candidates
    reco::Candidate* theGenZ = theGenAnalyzer->FindGenParticle(GenPVect, 23);
    reco::Candidate* theGenW = theGenAnalyzer->FindGenParticle(GenPVect, 24);
    // EWK corrections
    if(theGenZ) ZewkWeight = theGenAnalyzer->GetZewkWeight(theGenZ->pt());
    if(theGenW) WewkWeight = theGenAnalyzer->GetWewkWeight(theGenW->pt());
    
//    if(LheMap.find("lhePtZ")!=LheMap.end()) ZewkWeight = theGenAnalyzer->GetZewkWeight(LheMap["lhePtZ"]);
//    if(LheMap.find("lhePtW")!=LheMap.end()) WewkWeight = theGenAnalyzer->GetWewkWeight(LheMap["lhePtW"]);
    
    EventWeight *= ZewkWeight * WewkWeight;

    std::vector<int> LepIds = {12,14,16,-12,-14,-16};
    std::vector<int> HadIds = {1,2,3,4,5,-1,-2,-3,-4,-5};
    reco::GenParticle* theGenLep = theGenAnalyzer->FindGenParticleGenByIds(GenPVect, LepIds);
    reco::GenParticle* theGenHad = theGenAnalyzer->FindGenParticleGenByIds(GenPVect, HadIds);
    
    //Gen level plots and candidates
    //double GenZLepMass = -1.;
    //double GenZHadMass = -1.;
    //double GenZHadPt = -1.;
    //double GenXMass = -1.;
    double GenHadDR = -9.;
    double GenLepDR = -9.;
    //double GenZHadFatJetDR = -9.;
    bool isGenZZ = false;
    reco::Particle::LorentzVector p4GenZHad;
    if(theGenLep!=NULL && theGenHad!=NULL){
        const reco::Candidate* theGenZLep = theGenAnalyzer->FindMother(theGenLep);
        const reco::Candidate* theGenZHad = theGenAnalyzer->FindMother(theGenHad);
        for(unsigned int a = 0; a<=theGenZHad->numberOfDaughters(); a++) {
            if(theGenZHad!=NULL && theGenZHad->daughter(a)!=NULL && (theGenZHad->pdgId()==23 || theGenZHad->pdgId()==25) && (theGenZHad->daughter(a)->pdgId() == - theGenHad->pdgId())){
                GenHadDR = reco::deltaR(theGenHad->eta(),theGenHad->phi(),theGenZHad->daughter(a)->eta(),theGenZHad->daughter(a)->phi());
                break;
            }
        }
        for(unsigned int b = 0; b<=theGenZLep->numberOfDaughters(); b++) {
            if(theGenZLep!=NULL && theGenZLep->daughter(b)!=NULL && (theGenZLep->pdgId()==23 || theGenZLep->pdgId()==25) && (theGenZLep->daughter(b)->pdgId() == - theGenLep->pdgId())){
                GenLepDR = reco::deltaR(theGenLep->eta(),theGenLep->phi(),theGenZLep->daughter(b)->eta(),theGenZLep->daughter(b)->phi());
                break;
            }
        }
        if(theGenZLep!=NULL && theGenZLep->pdgId()==23 && theGenZHad!=NULL && (theGenZHad->pdgId()==23 || theGenZHad->pdgId()==25)) {
            Hist["g_ZLepMass"]->Fill(theGenZLep->mass(), EventWeight);
            Hist["g_ZLepPt"]->Fill(theGenZLep->pt(), EventWeight);
            Hist["g_LepPt"]->Fill(theGenLep->pt(), EventWeight);
            Hist["g_ZHadMass"]->Fill(theGenZHad->mass(), EventWeight);
            Hist["g_ZHadPt"]->Fill(theGenZHad->pt(), EventWeight);
            Hist["g_HadPt"]->Fill(theGenHad->pt(), EventWeight);
            Hist["g_HadEta"]->Fill(theGenHad->eta(), EventWeight);
            Hist["g_HadDR"]->Fill(GenHadDR, EventWeight);
            Hist["g_LepDR"]->Fill(GenLepDR, EventWeight);
            Hist["g_LepEta"]->Fill(theGenLep->eta(), EventWeight);
            Hist["g_ZZDR"]->Fill(reco::deltaR(theGenZHad->eta(),theGenZHad->phi(),theGenZLep->eta(),theGenZLep->phi()), EventWeight);
            Hist["g_ZZDPhi"]->Fill(reco::deltaPhi(theGenZHad->phi(),theGenZLep->phi()), EventWeight);
            Hist["g_LepHadDR"]->Fill(reco::deltaR(theGenHad->eta(),theGenHad->phi(),theGenLep->eta(),theGenLep->phi()), EventWeight);
            //GenZLepMass = theGenZLep->mass();
            //GenZHadMass = theGenZHad->mass();
            //GenZHadPt = theGenZHad->pt();
            isGenZZ = true;
            p4GenZHad = theGenZHad->p4();
        }
    }
    
    
    reco::Candidate* theGenX = theGenAnalyzer->FindGenParticleByIdAndStatus(GenPVect, 39, 62);
    if(!theGenX) theGenX = theGenAnalyzer->FindGenParticleByIdAndStatus(GenPVect, 36, 62);
    if(theGenX!=NULL && theGenLep!=NULL && theGenHad!=NULL && isGenZZ){
        Hist["g_XMass"]->Fill(theGenX->mass(), EventWeight);
        Hist["g_XMT"]->Fill(theGenX->mt(), EventWeight);
        Hist["g_XPt"]->Fill(theGenX->pt(), EventWeight);
        Hist["g_XRapidity"]->Fill(theGenX->rapidity(), EventWeight);
        //GenXMass = theGenX->mass();
    }
    

    // ---------- Trigger selections ----------
    // Dummy trigger
    //TriggerWeight*=theElectronAnalyzer->GetDoubleElectronTriggerSF(ElecVect.at(0), ElecVect.at(1));
    //TriggerWeight*=theMuonAnalyzer->GetDoubleMuonTriggerSF(MuonVect.at(0), MuonVect.at(1));
    
    
    
    Hist["a_nEvents"]->Fill(2., EventWeight);
    Hist["e_nEvents"]->Fill(2., EventWeight);
    Hist["m_nEvents"]->Fill(2., EventWeight);
    
    
    // Electron efficiency plots
    reco::GenParticle* genE1 = theGenAnalyzer->FindGenParticleGenByIds(GenPVect, std::vector<int>{-11});
    reco::GenParticle* genE2 = theGenAnalyzer->FindGenParticleGenByIds(GenPVect, std::vector<int>{+11});
    if(genE1!=NULL && genE2!=NULL) {
        int e1(-1), e2(-1);
        float dRll(-1.);
        dRll = reco::deltaR(genE1->eta(), genE1->phi(), genE2->eta(), genE2->phi());
        //pTll = (genE1->p4() + genE2->p4()).pt();
        
        for(unsigned int e = 0; e < ElecVect.size(); e++) {
            if(reco::deltaR(genE1->eta(), genE1->phi(), ElecVect[e].eta(), ElecVect[e].phi()) < 0.1 && fabs(1.-ElecVect[e].pt()/genE1->pt()) < 0.3) e1 = e;
            else if(((int)e)!=e1 && reco::deltaR(genE2->eta(), genE2->phi(), ElecVect[e].eta(), ElecVect[e].phi()) < 0.1 && fabs(1.-ElecVect[e].pt()/genE2->pt()) < 0.3) e2 = e;
        }
        if(e1 >= 0 && e2 >= 0) {
            Hist["e_nEvents"]->Fill(3., EventWeight);
            Hist["e_dR_reco"]->Fill(dRll);
            if(ElecVect[e1].pt() > 55. && ElecVect[e2].pt() > 20.) {
                Hist["e_nEvents"]->Fill(4., EventWeight);
                Hist["e_dR_pt"]->Fill(dRll);
                if(ElecVect[e1].charge() != ElecVect[e2].charge() && (ElecVect[e1].p4() + ElecVect[e2].p4()).mass() > 70 && (ElecVect[e1].p4() + ElecVect[e2].p4()).mass() < 110) {
                    Hist["e_nEvents"]->Fill(5., EventWeight);
                    Hist["e_dR_Z"]->Fill(dRll);
                    Hist["e_dR"]->Fill(dRll);
                    if(ElecVect[e1].userInt("isLoose") == 1 || ElecVect[e2].userInt("isLoose") == 1) {
                        Hist["e_nEvents"]->Fill(6., EventWeight);
                        Hist["e_dR_LooseId"]->Fill(dRll);
                        if(ElecVect[e1].userInt("isLoose") == 1 && ElecVect[e2].userInt("isLoose") == 1) {
                            Hist["e_nEvents"]->Fill(7., EventWeight);
                            Hist["e_dR_LooseLooseId"]->Fill(dRll);
                        }
                    }
                    if(ElecVect[e1].userInt("isVeto") == 1 && ElecVect[e2].userInt("isVeto") == 1) Hist["e_dR_VetoVetoId"]->Fill(dRll);
                    if(ElecVect[e1].userInt("isMedium") == 1 && ElecVect[e2].userInt("isMedium") == 1) Hist["e_dR_MediumMediumId"]->Fill(dRll);
                    if(ElecVect[e1].userInt("isTight") == 1 && ElecVect[e2].userInt("isTight") == 1) Hist["e_dR_TightTightId"]->Fill(dRll);
                }
            }
        }
    }
    
    
    // Muon efficiency plots
    reco::GenParticle* genM1 = theGenAnalyzer->FindGenParticleGenByIds(GenPVect, std::vector<int>{-13});
    reco::GenParticle* genM2 = theGenAnalyzer->FindGenParticleGenByIds(GenPVect, std::vector<int>{+13});
    if(genM1!=NULL && genM2!=NULL) {
        int m1(-1), m2(-1);
        float dRll(-1.);
        dRll = reco::deltaR(genM1->eta(), genM1->phi(), genM2->eta(), genM2->phi());
        //pTll = (genM1->p4() + genM2->p4()).pt();
        
        for(unsigned int m = 0; m < MuonVect.size(); m++) {
            if(reco::deltaR(genM1->eta(), genM1->phi(), MuonVect[m].eta(), MuonVect[m].phi()) < 0.1 && fabs(1.-MuonVect[m].pt()/genM1->pt()) < 0.3) m1 = m;
            else if(((int)m)!=m1 && reco::deltaR(genM2->eta(), genM2->phi(), MuonVect[m].eta(), MuonVect[m].phi()) < 0.1 && fabs(1.-MuonVect[m].pt()/genM2->pt()) < 0.3) m2 = m;
        }
        if(m1 >= 0 && m2 >= 0) {
            Hist["m_nEvents"]->Fill(3., EventWeight);
            Hist["m_dR_reco"]->Fill(dRll);
            if(MuonVect[m1].pt() > 55. && MuonVect[m2].pt() > 20.) {
                Hist["m_nEvents"]->Fill(4., EventWeight);
                Hist["m_dR_pt"]->Fill(dRll);
                if(MuonVect[m1].charge() != MuonVect[m2].charge() && (MuonVect[m1].p4() + MuonVect[m2].p4()).mass() > 70 && (MuonVect[m1].p4() + MuonVect[m2].p4()).mass() < 110) {
                    Hist["m_nEvents"]->Fill(5., EventWeight);
                    Hist["m_dR_Z"]->Fill(dRll);
                    Hist["m_dR"]->Fill(dRll);
                    if(MuonVect[m1].userInt("isHighPt") == 1 || MuonVect[m2].userInt("isHighPt") == 1) {
                        Hist["m_nEvents"]->Fill(6., EventWeight);
                        Hist["m_dR_HighptId"]->Fill(dRll);
                        if((MuonVect[m1].userInt("isHighPt")==1 && MuonVect[m2].userInt("isTrackerHighPt")==1) || (MuonVect[m2].userInt("isHighPt")==1 && MuonVect[m1].userInt("isTrackerHighPt")==1)) {
                            Hist["m_nEvents"]->Fill(7., EventWeight);
                            Hist["m_dR_HighptTrackerId"]->Fill(dRll);
                            if(MuonVect[m1].userFloat("trkIso") < 0.1 && MuonVect[m2].userFloat("trkIso") < 0.1) {
                                Hist["m_nEvents"]->Fill(8., EventWeight);
                                Hist["m_dR_HighptTrackerIdTrackerIso"]->Fill(dRll);
                            }
                        }
                    }
                    if(MuonVect[m1].userInt("isHighPt") == 1 && MuonVect[m2].userInt("isHighPt") == 1) Hist["m_dR_HighptHighptId"]->Fill(dRll);
                    if(MuonVect[m1].userInt("isTight") == 1 && MuonVect[m2].userInt("isTight") == 1) Hist["m_dR_TightTightId"]->Fill(dRll);
                    if(MuonVect[m1].userInt("isMedium") == 1 && MuonVect[m2].userInt("isMedium") == 1) Hist["m_dR_MediumMediumId"]->Fill(dRll);
                    if(MuonVect[m1].userInt("isLoose") == 1 && MuonVect[m2].userInt("isLoose") == 1) Hist["m_dR_LooseLooseId"]->Fill(dRll);
                }
            }
        }
    }
  
    
    // ---------- Event Variables ----------
    
    // Max b-tagged jet in the event
    for(unsigned int i = 2; i < JetsVect.size(); i++) if(JetsVect[i].bDiscriminator(JetPSet.getParameter<std::string>("btag")) > MaxJetBTag) MaxJetBTag = JetsVect[i].bDiscriminator(JetPSet.getParameter<std::string>("btag"));
    // Max b-tagged jet in the event
    for(unsigned int i = 0; i < JetsVect.size(); i++) if(FatJetsVect.size() > 0 && JetsVect[i].bDiscriminator(JetPSet.getParameter<std::string>("btag")) > MaxFatJetBTag && deltaR(FatJetsVect.at(0), JetsVect[i])>0.8) MaxFatJetBTag = JetsVect[i].bDiscriminator(JetPSet.getParameter<std::string>("btag"));
    
    for(unsigned int i = 0; i < JetsVect.size(); i++) if(fabs(reco::deltaPhi(JetsVect[i].phi(), MET.phi())) < MinJetMetDPhi) MinJetMetDPhi = fabs(reco::deltaPhi(JetsVect[i].phi(), MET.phi()));
    
    // Jet variables
    theJetAnalyzer->AddVariables(JetsVect, MET);
    theFatJetAnalyzer->AddVariables(FatJetsVect, MET);
    // Leptons
    theElectronAnalyzer->AddVariables(ElecVect, MET);
    theMuonAnalyzer->AddVariables(MuonVect, MET);
    
    // Highest CSV bjet in the event

    // ---------- Print Summary ----------
    if(Verbose) {
        std::cout << " --- Event n. " << iEvent.id().event() << ", lumi " << iEvent.luminosityBlock() << ", run " << iEvent.id().run() << ", weight " << EventWeight << std::endl;
        std::cout << "number of electrons: " << ElecVect.size() << std::endl;
        for(unsigned int i = 0; i < ElecVect.size(); i++) std::cout << "  electron [" << i << "]\tpt: " << ElecVect[i].pt() << "\teta: " << ElecVect[i].eta() << "\tphi: " << ElecVect[i].phi() << "\tmass: " << ElecVect[i].mass() << "\tcharge: " << ElecVect[i].charge() << std::endl;
        std::cout << "number of muons:     " << MuonVect.size() << std::endl;
        for(unsigned int i = 0; i < MuonVect.size(); i++) std::cout << "  muon     [" << i << "]\tpt: " << MuonVect[i].pt() << "\teta: " << MuonVect[i].eta() << "\tphi: " << MuonVect[i].phi() << "\tmass: " << MuonVect[i].mass() << "\tcharge: " << MuonVect[i].charge() << std::endl;
        std::cout << "number of AK4 jets:  " << JetsVect.size() << std::endl;
        for(unsigned int i = 0; i < JetsVect.size(); i++) std::cout << "  AK4 jet  [" << i << "]\tpt: " << JetsVect[i].pt() << "\teta: " << JetsVect[i].eta() << "\tphi: " << JetsVect[i].phi() << "\tmass: " << JetsVect[i].mass() << std::endl;
        std::cout << "number of AK8 jets:  " << FatJetsVect.size() << std::endl;
        for(unsigned int i = 0; i < FatJetsVect.size(); i++) std::cout << "  AK8 jet  [" << i << "]\tpt: " << FatJetsVect[i].pt() << "\teta: " << FatJetsVect[i].eta() << "\tphi: " << FatJetsVect[i].phi() << "\tmass: " << FatJetsVect[i].mass() << std::endl;
        std::cout << "Missing energy:      " << MET.pt() << std::endl;
    }

    
    // ---------- Fill objects ----------
    if(Verbose) std::cout << " - Filling objects" << std::endl;
    if(isZtoEE || isWtoEN) for(unsigned int i = 0; i < Leptons.size() && i < ElecVect.size(); i++) ObjectsFormat::FillElectronType(Leptons[i], &ElecVect[i], isMC);
    else if(isZtoMM || isWtoMN) for(unsigned int i = 0; i < Leptons.size() && i < MuonVect.size(); i++) ObjectsFormat::FillMuonType(Leptons[i], &MuonVect[i], isMC);
    else if(isTtoEM && Leptons.size() >= 2) {
        if(ElecVect[0].pt() > MuonVect[0].pt()) {
            ObjectsFormat::FillElectronType(Leptons[0], &ElecVect[0], isMC);
            ObjectsFormat::FillMuonType(Leptons[1], &MuonVect[0], isMC);
        }
        else {
            ObjectsFormat::FillMuonType(Leptons[0], &MuonVect[0], isMC);
            ObjectsFormat::FillElectronType(Leptons[1], &ElecVect[0], isMC);
        }
    }
    // clear previous vector
    Jets.clear();
    for(unsigned int i = 0; i < JetsVect.size(); i++) {
      Jets.emplace_back();
      ObjectsFormat::FillJetType(Jets[i], &JetsVect[i], isMC);
    }
    ObjectsFormat::FillMEtType(MEt, &MET, isMC);
    
    // Fill tree
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
void HHAnalyzer::beginJob() {
    
    // Object objects are created only one in the begin job. The reference passed to the branch has to be the same
    for(int i = 0; i < WriteNElectrons; i++) Electrons.push_back( LeptonType() );
    for(int i = 0; i < WriteNMuons; i++) Muons.push_back( LeptonType() );
    for(int i = 0; i < WriteNLeptons; i++) Leptons.push_back( LeptonType() );
    for(int i = 0; i < WriteNFatJets; i++) FatJets.push_back( FatJetType() );
    
    // Create Tree and set Branches
    tree=fs->make<TTree>("tree", "tree");
    tree->Branch("isMC", &isMC, "isMC/O");
    tree->Branch("EventNumber", &EventNumber, "EventNumber/L");
    tree->Branch("LumiNumber", &LumiNumber, "LumiNumber/L");
    tree->Branch("RunNumber", &RunNumber, "RunNumber/L");
    tree->Branch("EventWeight", &EventWeight, "EventWeight/F");
    tree->Branch("FacWeightUp", &FacWeightUp, "FacWeightUp/F");
    tree->Branch("FacWeightDown", &FacWeightDown, "FacWeightDown/F");
    tree->Branch("RenWeightUp", &RenWeightUp, "RenWeightUp/F");
    tree->Branch("RenWeightDown", &RenWeightDown, "RenWeightDown/F");
    tree->Branch("ScaleWeightUp", &ScaleWeightUp, "ScaleWeightUp/F");
    tree->Branch("ScaleWeightDown", &ScaleWeightDown, "ScaleWeightDown/F");
    tree->Branch("StitchWeight", &StitchWeight, "StitchWeight/F");
    tree->Branch("ZewkWeight", &ZewkWeight, "ZewkWeight/F");
    tree->Branch("WewkWeight", &WewkWeight, "WewkWeight/F");
    tree->Branch("PUWeight", &PUWeight, "PUWeight/F");
    tree->Branch("TriggerWeight", &TriggerWeight, "TriggerWeight/F");
    tree->Branch("LeptonWeight", &LeptonWeight, "LeptonWeight/F");
    
    // Set trigger branches
    for(auto it = TriggerMap.begin(); it != TriggerMap.end(); it++) tree->Branch(it->first.c_str(), &(it->second), (it->first+"/O").c_str());
    
    // Analysis variables
    tree->Branch("isZtoEE", &isZtoEE, "isZtoEE/O");
    tree->Branch("isZtoMM", &isZtoMM, "isZtoMM/O");
    tree->Branch("isTtoEM", &isTtoEM, "isTtoEM/O");
    tree->Branch("isWtoEN", &isWtoEN, "isWtoEN/O");
    tree->Branch("isWtoMN", &isWtoMN, "isWtoMN/O");
    tree->Branch("isZtoNN", &isZtoNN, "isZtoNN/O");
    tree->Branch("isMerged", &isMerged, "isMerged/O");
    tree->Branch("isResolved", &isResolved, "isResolved/O");
    
    // Objects
    tree->Branch("nPV", &nPV, "nPV/I");
    tree->Branch("nElectrons", &nElectrons, "nElectrons/I");
    tree->Branch("nVetoElectrons", &nVetoElectrons, "nVetoElectrons/I");
    tree->Branch("nMuons", &nMuons, "nMuons/I");
    tree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
    tree->Branch("nJets", &nJets, "nJets/I");
    tree->Branch("nFatJets", &nFatJets, "nFatJets/I");
    tree->Branch("nBTagJets", &nBTagJets, "nBTagJets/I");
    
    tree->Branch("MaxJetBTag", &MaxJetBTag, "MaxJetBTag/F");
    tree->Branch("MaxFatJetBTag", &MaxFatJetBTag, "MaxFatJetBTag/F");
    tree->Branch("MinJetMetDPhi", &MinJetMetDPhi, "MinJetMetDPhi/F");
    tree->Branch("Chi2", &Chi2, "Chi2/F");
    // Angular variables
    tree->Branch("CosThetaStar", &CosThetaStar, "CosThetaStar/F");
    tree->Branch("CosTheta1", &CosTheta1, "CosTheta1/F");
    tree->Branch("CosTheta2", &CosTheta2, "CosTheta2/F");
    tree->Branch("Phi", &Phi, "Phi/F");
    tree->Branch("Phi1", &Phi1, "Phi1/F");
    // Mass recoil formula
    tree->Branch("massRecoilFormula", &massRecoilFormula, "massRecoilFormula/F");
  
    // Set Branches for objects
    for(int i = 0; i < WriteNElectrons; i++) tree->Branch(("Electron"+std::to_string(i+1)).c_str(), &(Electrons[i].pt), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNMuons; i++) tree->Branch(("Muon"+std::to_string(i+1)).c_str(), &(Muons[i].pt), ObjectsFormat::ListLeptonType().c_str());
    for(int i = 0; i < WriteNLeptons; i++) tree->Branch(("Lepton"+std::to_string(i+1)).c_str(), &(Leptons[i].pt), ObjectsFormat::ListLeptonType().c_str());
    // save vector of jets
    tree->Branch("Jets", &(Jets), 64000, 2);
    for(int i = 0; i < WriteNFatJets; i++) tree->Branch(("FatJet"+std::to_string(i+1)).c_str(), &(FatJets[i].pt), ObjectsFormat::ListFatJetType().c_str());
    tree->Branch("MEt", &MEt.pt, ObjectsFormat::ListMEtType().c_str());
    tree->Branch("V", &V.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("X", &X.pt, ObjectsFormat::ListCandidateType().c_str());
    /*tree->Branch("H", &H.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("A", &A.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("kH", &kH.pt, ObjectsFormat::ListLorentzType().c_str());
    tree->Branch("kA", &kA.pt, ObjectsFormat::ListLorentzType().c_str());
    tree->Branch("HMerged", &HMerged.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("XMerged", &XMerged.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("HResolved", &HResolved.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("HResolvedPt", &HResolvedPt.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("HResolvedHpt", &HResolvedHpt.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("HResolvedDZ", &HResolvedDZ.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("HResolvedDR", &HResolvedDR.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("XResolved", &XResolved.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("XResolvedPt", &XResolvedPt.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("XResolvedHpt", &XResolvedHpt.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("XResolvedDZ", &XResolvedDZ.pt, ObjectsFormat::ListCandidateType().c_str());
    tree->Branch("XResolvedDR", &XResolvedDR.pt, ObjectsFormat::ListCandidateType().c_str());
    */
    // -------------------    
    /*
    // Create Tree for alpha and set Branches    
    treealpha=fs->make<TTree>("treealpha", "treealpha");

    treealpha->Branch("isMC", &isMC, "isMC/O");
    treealpha->Branch("EventWeight", &EventWeight, "EventWeight/F");
    
    treealpha->Branch("EventNumber", &EventNumber, "EventNumber/L");
    treealpha->Branch("LumiNumber", &LumiNumber, "LumiNumber/L");
    treealpha->Branch("RunNumber", &RunNumber, "RunNumber/L");

    // Set trigger branches
    for(auto it = TriggerMap.begin(); it != TriggerMap.end(); it++) treealpha->Branch(it->first.c_str(), &(it->second), (it->first+"/O").c_str());
    
    // Analysis variables
    treealpha->Branch("isZtoEE", &isZtoEE, "isZtoEE/O");
    treealpha->Branch("isZtoMM", &isZtoMM, "isZtoMM/O");
    
    // Lepton1
    treealpha->Branch("Lepton1_isMuon", &Lepton1_isMuon, "Lepton1_isMuon/O");
    treealpha->Branch("Lepton1_isElectron", &Lepton1_isElectron, "Lepton1_isElectron/O");
    treealpha->Branch("Lepton1_isLoose", &Lepton1_isLoose, "Lepton1_isLoose/O");
    treealpha->Branch("Lepton1_isHighPt", &Lepton1_isHighPt, "Lepton1_isHighPt/O");
    treealpha->Branch("Lepton1_isTrackerHighPt", &Lepton1_isTrackerHighPt, "Lepton1_isTrackerHighPt/O");
    treealpha->Branch("Lepton1_pt", &Lepton1_pt, "Lepton1_pt/F");
    treealpha->Branch("Lepton1_trkIso", &Lepton1_trkIso, "Lepton1_trkIso/F");

    // Lepton2
    treealpha->Branch("Lepton2_isMuon", &Lepton2_isMuon, "Lepton2_isMuon/O");
    treealpha->Branch("Lepton2_isElectron", &Lepton2_isElectron, "Lepton2_isElectron/O");
    treealpha->Branch("Lepton2_isLoose", &Lepton2_isLoose, "Lepton2_isLoose/O");
    treealpha->Branch("Lepton2_isHighPt", &Lepton2_isHighPt, "Lepton2_isHighPt/O");
    treealpha->Branch("Lepton2_isTrackerHighPt", &Lepton2_isTrackerHighPt, "Lepton2_isTrackerHighPt/O");
    treealpha->Branch("Lepton2_pt", &Lepton2_pt, "Lepton2_pt/F");
    treealpha->Branch("Lepton2_trkIso", &Lepton2_trkIso, "Lepton2_trkIso/F");

    // MET        
    treealpha->Branch("MEt_pt", &MEt_pt, "MEt_pt/F");

    // V        
    treealpha->Branch("V_pt", &V_pt, "V_pt/F");
    treealpha->Branch("V_mass", &V_mass, "V_mass/F");
    treealpha->Branch("V_tmass", &V_tmass, "V_tmass/F");

    // X        
    treealpha->Branch("X_pt", &X_pt, "X_pt/F");
    treealpha->Branch("X_dPhi", &X_dPhi, "X_dPhi/F");
    treealpha->Branch("X_mass", &X_mass, "X_mass/F");
    treealpha->Branch("X_tmass", &X_tmass, "X_tmass/F");

    // FatJet1
    treealpha->Branch("FatJet1_pt", &FatJet1_pt, "FatJet1_pt/F");
    treealpha->Branch("FatJet1_prunedMass", &FatJet1_prunedMass, "FatJet1_prunedMass/F");
    treealpha->Branch("FatJet1_softdropMass", &FatJet1_softdropMass, "FatJet1_softdropMass/F");
    treealpha->Branch("FatJet1_softdropPuppiMass", &FatJet1_softdropPuppiMass, "FatJet1_softdropPuppiMass/F");
    treealpha->Branch("FatJet1_prunedMassCorr", &FatJet1_prunedMassCorr, "FatJet1_prunedMassCorr/F");
    treealpha->Branch("FatJet1_softdropMassCorr", &FatJet1_softdropMassCorr, "FatJet1_softdropMassCorr/F");
    treealpha->Branch("FatJet1_softdropPuppiMassCorr", &FatJet1_softdropPuppiMassCorr, "FatJet1_softdropPuppiMassCorr/F");
    treealpha->Branch("FatJet1_chsTau21", &FatJet1_chsTau21, "FatJet1_chsTau21/F");
    treealpha->Branch("FatJet1_puppiTau21", &FatJet1_puppiTau21, "FatJet1_puppiTau21/F");
    treealpha->Branch("FatJet1_ddtTau21", &FatJet1_ddtTau21, "FatJet1_ddtTau21/F");
    */
}

// ------------ method called once each job just after ending the event loop  ------------
void HHAnalyzer::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HHAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(HHAnalyzer);
