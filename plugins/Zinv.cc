// -*- C++ -*-
//
// Package:    Analysis/Zinv
// Class:      Zinv
// 
/**\class Zinv Zinv.cc Analysis/Zinv/plugins/Zinv.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alberto Zucchetta
//         Created:  Thu, 28 Apr 2016 08:28:54 GMT
//
//

#include "Zinv.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Zinv::Zinv(const edm::ParameterSet& iConfig):
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
    TFileDirectory tauDir=fs->mkdir("Taus/");
    TFileDirectory phoDir=fs->mkdir("Photons/");
    TFileDirectory jetDir=fs->mkdir("Jets/");
    TFileDirectory kinDir=fs->mkdir("Kin/");
    
    // Make TH1F
    //std::vector<std::string> nLabels={"All", "Trigger", "Iso Lep #geq 2", "Z cand ", "Jets #geq 2", "Z mass ", "h mass ", "Top veto", "bJets #geq 1", "bJets #geq 2"};
    std::vector<std::string> nLabels={"All", "METcut", "FatJetLoose", "noEle", "noMu", "noTau", "noPhotons", "Trigger"};
    
    int nbins;
    float min, max;
    std::string name, title, opt;
    
    ifstream histFile(HistFile);
    if(!histFile.is_open()) {
        throw cms::Exception("Zinv Analyzer", HistFile + " file not found");
    }
    while(histFile >> name >> title >> nbins >> min >> max >> opt) {
        if(name.find('#')==std::string::npos) {
            while(title.find("~")!=std::string::npos) title=title.replace(title.find("~"), 1, " "); // Remove ~
            if(name.substr(0, 2)=="a_") Hist[name] = allDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max); //.substr(2)
            if(name.substr(0, 2)=="g_") Hist[name] = genDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            if(name.substr(0, 2)=="e_") Hist[name] = eleDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            if(name.substr(0, 2)=="m_") Hist[name] = muoDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            if(name.substr(0, 2)=="t_") Hist[name] = tauDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            if(name.substr(0, 2)=="p_") Hist[name] = phoDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            if(name.substr(0, 2)=="j_") Hist[name] = jetDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            if(name.substr(0, 2)=="k_") Hist[name] = kinDir.make<TH1F>(name.c_str(), title.c_str(), nbins, min, max);
            Hist[name]->Sumw2();
            Hist[name]->SetOption(opt.c_str());
            // Particular histograms
            if(name=="a_nEvents" || name=="e_nEvents" || name=="m_nEvents") for(unsigned int i=0; i<nLabels.size(); i++) Hist[name]->GetXaxis()->SetBinLabel(i+1, nLabels[i].c_str());
        }
    }
    histFile.close();

    //TH2F booking hardcoded
    MetDPhiNoMuPreTrig_MuSel = fs->make<TH2F>("MetDPhiNoMuPreTrig_MuSel","MetDPhiNoMuPreTrig_MuSel",500,0,5000,100,-5,5);
    MetDPhiWithMuPreTrig_MuSel = fs->make<TH2F>("MetDPhiWithMuPreTrig_MuSel","MetDPhiWithMuPreTrig_MuSel",500,0,5000,100,-5,5);
    MetDPhiNoMuPostTrig_MuSel = fs->make<TH2F>("MetDPhiNoMuPostTrig_MuSel","MetDPhiNoMuPostTrig_MuSel",500,0,5000,100,-5,5);
    MetDPhiWithMuPostTrig_MuSel = fs->make<TH2F>("MetDPhiWithMuPostTrig_MuSel","MetDPhiWithMuPostTrig_MuSel",500,0,5000,100,-5,5);
    //MetDPhiPreTrig_MuSel = fs->make<TH2F>("MetDPhiPreTrig_MuSel","MetDPhiPreTrig_MuSel",500,0,5000,100,-5,5);
    //MetDPhiPostTrig_MuSel = fs->make<TH2F>("MetDPhiPostTrig_MuSel","MetDPhiPostTrig_MuSel",500,0,5000,100,-5,5);

    MetDPhiPreTrig_EleSel = fs->make<TH2F>("MetDPhiPreTrig_EleSel","MetDPhiPreTrig_EleSel",500,0,5000,100,-5,5);
    MetDPhiPostTrig_EleSel = fs->make<TH2F>("MetDPhiPostTrig_EleSel","MetDPhiPostTrig_EleSel",500,0,5000,100,-5,5);
    
    std::cout << "---------- STARTING ----------" << std::endl;
}


Zinv::~Zinv() {
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
void Zinv::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    isMC = !iEvent.isRealData();
    EventNumber = iEvent.id().event();
    LumiNumber = iEvent.luminosityBlock();
    RunNumber = iEvent.id().run();
    
    EventWeight = StitchWeight = ZewkWeight = WewkWeight = PUWeight = TriggerWeight = LeptonWeight = 1.;
    FacWeightUp = FacWeightDown = RenWeightUp = RenWeightDown = ScaleWeightUp = ScaleWeightDown = 1.;
    isZtoEE = isZtoMM = isTtoEM = isWtoEN = isWtoMN = isZtoNN = false;
    nPV = nElectrons = nMuons = nTaus = nPhotons = nJets = nFatJets = nBTagJets = 1;
    nVetoElectrons = nLooseElectrons = nLooseMuons = nTrackerHighPtMuons = 0;
    MaxJetBTag = MaxFatJetBTag = Chi2 = -1.;
    MinJetMetDPhi = 10.;
    massRecoilFormula = -1.;
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
    ObjectsFormat::ResetCandidateType(X);
    
    Hist["a_nEvents"]->Fill(1., EventWeight);
    
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
    //EventWeight *= TriggerWeight;
    
    // Electrons
    std::vector<pat::Electron> ElecVect = theElectronAnalyzer->FillElectronVector(iEvent);
    nElectrons = ElecVect.size();
    for(unsigned int i=0; i<ElecVect.size(); i++){
        if(ElecVect.at(i).userInt("isVeto")==1) nVetoElectrons++;
        if(ElecVect.at(i).userInt("isLoose")==1) nLooseElectrons++;
    }
    // Muons
    std::vector<pat::Muon> MuonVect = theMuonAnalyzer->FillMuonVector(iEvent);
    nMuons = MuonVect.size();
    float nTightMuons = 0.;
    for(unsigned int i=0; i<MuonVect.size(); i++){
        if(MuonVect.at(i).userInt("isLoose")==1) nLooseMuons++;
        if(MuonVect.at(i).userInt("isTight")==1) nTightMuons++;
        if(MuonVect.at(i).userInt("isTrackerHighPt")==1) nTrackerHighPtMuons++;
    }
    // Taus
    std::vector<pat::Tau> TauVect = theTauAnalyzer->FillTauVector(iEvent);
    theTauAnalyzer->CleanTausFromMuons(TauVect, MuonVect, 0.4);
    theTauAnalyzer->CleanTausFromElectrons(TauVect, ElecVect, 0.4);
    nTaus = TauVect.size();
    // Photons
    std::vector<pat::Photon> PhotonVect = thePhotonAnalyzer->FillPhotonVector(iEvent);
    if(TriggerMap.find("HLT_DoublePhoton60_v") != TriggerMap.end() && TriggerMap["HLT_DoublePhoton60_v"]) thePhotonAnalyzer->PlotPhotons(PhotonVect, Hist, EventWeight);
    nPhotons = PhotonVect.size();
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
    //MET trigger efficiency for muons ad electrons
    float metNoMupx = MET.px();
    float metNoMupy = MET.py();
    bool METtrigNoMu = false;
    bool METtrig = false;
    METtrigNoMu = TriggerMap["HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v"] || TriggerMap["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"];
    //METtrig = TriggerMap["HLT_PFMET170_NotCleaned_v"] || TriggerMap["HLT_PFMET170_NoiseCleaned_v"] || TriggerMap["HLT_PFMET170_JetIdCleaned_v"] || TriggerMap["HLT_PFMET170_HBHECleaned_v"] || TriggerMap["HLT_PFMET170_BeamHaloCleaned_v"];
    METtrig = TriggerMap["HLT_PFMET170_NoiseCleaned_v"] || TriggerMap["HLT_PFMET170_JetIdCleaned_v"] || TriggerMap["HLT_PFMET170_HBHECleaned_v"];

    if(nFatJets>0){
        if(nMuons>0){
            if(TriggerMap["HLT_Mu45_eta2p1_v"] && MuonVect.at(0).userInt("isTight") && MuonVect.at(0).userFloat("pfIso04")<0.15 && MuonVect.at(0).pt()>55){
                for(int m=0; m<nMuons; m++){
	            metNoMupx += MuonVect.at(m).px();
	            metNoMupy += MuonVect.at(m).py();
		}
	        reco::Particle::LorentzVector metNoMup4(metNoMupx, metNoMupy, 0, sqrt(pow(metNoMupx,2) + pow(metNoMupy,2)) );
	        MET.addUserFloat("ptNoMu",metNoMup4.pt());
	        MET.addUserFloat("phiNoMu",metNoMup4.phi());
	        Hist["a_mu_den_trig"]->Fill(MET.pt(), EventWeight);
	        Hist["a_mu_den_trigNoMu"]->Fill(MET.userFloat("ptNoMu"), EventWeight);
		MetDPhiNoMuPreTrig_MuSel->Fill(MET.userFloat("ptNoMu"), reco::deltaPhi(MuonVect.at(0).phi(),MET.userFloat("phiNoMu")), EventWeight);
		MetDPhiWithMuPreTrig_MuSel->Fill(MET.pt(), reco::deltaPhi(MuonVect.at(0).phi(),MET.phi()), EventWeight);
	        if(METtrigNoMu){
		    Hist["a_mu_num_trigNoMu"]->Fill(MET.userFloat("ptNoMu"), EventWeight);
		    MetDPhiNoMuPostTrig_MuSel->Fill(MET.userFloat("ptNoMu"), reco::deltaPhi(MuonVect.at(0).phi(),MET.userFloat("phiNoMu")), EventWeight);
		}
	        if(METtrig){
		    Hist["a_mu_num_trig"]->Fill(MET.pt(), EventWeight);
		    MetDPhiWithMuPostTrig_MuSel->Fill(MET.pt(), reco::deltaPhi(MuonVect.at(0).phi(),MET.phi()), EventWeight);
		}
            }
        }
        if(nElectrons>0){
            if(TriggerMap["HLT_Ele27_WPLoose_Gsf_v"] && ElecVect.at(0).userInt("isTight") && ElecVect.at(0).pt()>55){
	        Hist["a_ele_den_trig"]->Fill(MET.pt(), EventWeight);
		MetDPhiPreTrig_EleSel->Fill(MET.pt(), reco::deltaPhi(ElecVect.at(0).phi(),MET.phi()), EventWeight);
	        if(METtrig || METtrigNoMu) {
		    Hist["a_ele_num_trig"]->Fill(MET.pt(), EventWeight);
		    MetDPhiPostTrig_EleSel->Fill(MET.pt(), reco::deltaPhi(ElecVect.at(0).phi(),MET.phi()), EventWeight);
		}
            }
        }
    }
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
    
    
    
    
    // -----------------------------------
    //           VECTOR BOSON
    // -----------------------------------
    
    // Categorization depending on the number of leptons

    // ----------- Z TO NEUTRINOS ---------------
    if(MET.pt() > 100.) {
        isZtoNN=true;
        if(Verbose) std::cout << " - No charged leptons" << std::endl;
    }

    // ---------- Reconstruct V candidate ----------
    
    pat::CompositeCandidate theV;
    if(isZtoNN) {
      if(Verbose) std::cout << " - Try to reconstruct Z -> nn" << std::endl;
      Hist["a_nEvents"]->Fill(2., EventWeight);
      isZtoNN=true;
      theV.addDaughter(MET);
      addP4.set(theV);
      if(FatJetsVect.size()>0 && FatJetsVect.at(0).hasUserInt("isLoose")){
	Hist["a_nEvents"]->Fill(3., EventWeight);
	if(nVetoElectrons==0) {
	  Hist["a_nEvents"]->Fill(4., EventWeight);
	  if(nLooseMuons==0) {
	    Hist["a_nEvents"]->Fill(5., EventWeight);
	    if(nTaus==0) {
	      Hist["a_nEvents"]->Fill(6., EventWeight);
	      if(nPhotons==0){
		Hist["a_nEvents"]->Fill(7., EventWeight);
		if(METtrig || METtrigNoMu){
		  Hist["a_nEvents"]->Fill(8., EventWeight);
		}
	      }
	    }
	  }
	}
      }
	    //if(Verbose) std::cout << " - Possible Z inv candidate" << std::endl;    
	    /*
            if(isMC) {
                /// FIXME -> APPLYING THE SF FOR Met trigger HARDCODED <- FIXME ///
                TriggerWeight *= theJetAnalyzer->GetMetTriggerPFMETNoMu90OrPFMETNoMu120OrPFMET170SF(MET);
            }
	    */
    }
    
    if(not(isZtoNN) && not(isZtoEE) && not(isZtoMM)) { if(Verbose) std::cout << " - No reconstructible V candidate" << std::endl; return; }
    // Update event weight with lepton selections
    EventWeight *= LeptonWeight;
    // Update event weight with met trigger weight
    //EventWeight *= TriggerWeight;
    
    

    // -----------------------------------
    //           HADRONIC BOSON
    // -----------------------------------
    
    if(Verbose) std::cout << " - Cleaning FatJet" << std::endl;
    for(unsigned int j = 0; j < FatJetsVect.size(); ) {
        if((isZtoEE || isZtoMM) && (deltaR(FatJetsVect[j], *theV.daughter(0)) < 0.8 || deltaR(FatJetsVect[j], *theV.daughter(1)) < 0.8)) FatJetsVect.erase(FatJetsVect.begin() + j);
        else j++;
    }
    if(FatJetsVect.size() < 1) {if(Verbose) std::cout << " - N fat jets < 1" << std::endl; return;}

    // ------------------------------
    // ------- RES CANDIDATE --------
    // ------------------------------
    pat::CompositeCandidate theX;
    theX.addDaughter(theV);
    theX.addDaughter(FatJetsVect.at(0));
    addP4.set(theX);
    if(isZtoNN && theX.numberOfDaughters()>0) {
        float D = pow(FatJetsVect.at(0).energy(),2) - (pow(FatJetsVect.at(0).mass(),2) - pow(90.18,2));
        if(D>0) {
            float s1 = FatJetsVect.at(0).energy() + sqrt(D);
            float s2 = FatJetsVect.at(0).energy() - sqrt(D);
            massRecoilFormula = (fabs(s1)>fabs(s2))? s1 : s2;
        }
        else massRecoilFormula = FatJetsVect.at(0).energy();
    }
    if(Verbose) std::cout << " - Candidate built" << std::endl;
    
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
        std::cout << "isZtoEE: " << isZtoEE << ", is ZtoMM: " << isZtoMM << ", isZtoNN: " << isZtoNN << std::endl;
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
        std::cout << "V leptonic mass:     " << theV.mass() << std::endl;
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
    for(unsigned int i = 0; i < Taus.size() && i < TauVect.size(); i++) ObjectsFormat::FillTauType(Taus[i], &TauVect[i], isMC);
    for(unsigned int i = 0; i < Photons.size() && i < PhotonVect.size(); i++) ObjectsFormat::FillPhotonType(Photons[i], &PhotonVect[i], isMC);
    for(unsigned int i = 0; i < Jets.size() && i < JetsVect.size(); i++) ObjectsFormat::FillJetType(Jets[i], &JetsVect[i], isMC);
    for(unsigned int i = 0; i < FatJets.size() && i < FatJetsVect.size(); i++) ObjectsFormat::FillFatJetType(FatJets[i], &FatJetsVect[i], isMC);
    ObjectsFormat::FillMEtType(MEt, &MET, isMC);
    ObjectsFormat::FillCandidateType(V, &theV, isMC);
    ObjectsFormat::FillCandidateType(X, &theX, isMC);
    if(V.pt < 100.) return;
    
    // Fill tree
    tree->Fill();
    
}


// ------------ method called once each job just before starting event loop  ------------
void Zinv::beginJob() {
    
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
    tree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
    tree->Branch("nMuons", &nMuons, "nMuons/I");
    tree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
    tree->Branch("nTrackerHighPtMuons", &nTrackerHighPtMuons, "nTrackerHighPtMuons/I");
    tree->Branch("nTaus", &nTaus, "nTaus/I");
    tree->Branch("nPhotons", &nPhotons, "nPhotons/I");
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
    for(int i = 0; i < WriteNTaus; i++) tree->Branch(("Tau"+std::to_string(i+1)).c_str(), &(Taus[i].pt), ObjectsFormat::ListTauType().c_str());
    for(int i = 0; i < WriteNPhotons; i++) tree->Branch(("Photon"+std::to_string(i+1)).c_str(), &(Photons[i].pt), ObjectsFormat::ListPhotonType().c_str());
    for(int i = 0; i < WriteNJets; i++) tree->Branch(("Jet"+std::to_string(i+1)).c_str(), &(Jets[i].pt), ObjectsFormat::ListJetType().c_str());
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
void Zinv::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Zinv::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(Zinv);
