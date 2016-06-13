//#include "RecoilCorrector.hh" // From: https://github.com/cms-met/MetTools/tree/master/RecoilCorrections

#include "JetAnalyzer.h"


JetAnalyzer::JetAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    JetToken(CColl.consumes<std::vector<pat::Jet> >(PSet.getParameter<edm::InputTag>("jets"))),
    MetToken(CColl.consumes<std::vector<pat::MET> >(PSet.getParameter<edm::InputTag>("met"))),
    CorToken(CColl.consumes<reco::JetCorrector>(PSet.getParameter<edm::InputTag>("corrector"))),
    QGToken(CColl.consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"))),
    JetId(PSet.getParameter<int>("jetid")),
    Jet1Pt(PSet.getParameter<double>("jet1pt")),
    Jet2Pt(PSet.getParameter<double>("jet2pt")),
    JetEta(PSet.getParameter<double>("jeteta")),
    AddQG(PSet.getParameter<bool>("addQGdiscriminator")),
    RecalibrateJets(PSet.getParameter<bool>("recalibrateJets")),
    RecalibrateMass(PSet.getParameter<bool>("recalibrateMass")),
    UseReshape(PSet.getParameter<bool>("reshapeBTag")),
    BTag(PSet.getParameter<std::string>("btag")),
    Jet1BTag(PSet.getParameter<int>("jet1btag")),
    Jet2BTag(PSet.getParameter<int>("jet2btag")),
    BTagDB(PSet.getParameter<std::string>("btagDB")),
    UseRecoil(PSet.getParameter<bool>("metRecoil")),
    RecoilMCFile(PSet.getParameter<std::string>("metRecoilMC")),
    RecoilDataFile(PSet.getParameter<std::string>("metRecoilData"))
{
  
    isJESFile=false;

//    JESFile=new TFile("data/JESUncertainty.root", "READ");
//    JESFile->cd();
//    if(!JESFile->IsZombie()) {
//      hist=(TH2F*)JESFile->Get("test/AK5PFchs");
//      isJESFile=true;
//    }
    
    
    // BTag calibrator
    if(UseReshape) {
        calib           = new BTagCalibration("CSVv2", BTagDB);
        
        reader          = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","central");
        reader_up_jes   = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","up_jes");
        reader_down_jes = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","down_jes");
    //    reader_up_lf        = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","up_lf");
    //    reader_up_hfstats1  = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","up_hfstats1");
    //    reader_up_hfstats2  = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","up_hfstats2");
    //    reader_up_cferr1    = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","up_cferr1");
    //    reader_up_cferr2    = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","up_cferr2");
    //    reader_down_lf        = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","down_lf");
    //    reader_down_hfstats1  = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","down_hfstats1");
    //    reader_down_hfstats2  = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","down_hfstats2");
    //    reader_down_cferr1    = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","down_cferr1");
    //    reader_down_cferr2    = new BTagCalibrationReader(calib,BTagEntry::OP_RESHAPING,"iterativefit","down_cferr2");
    }
    // Recoil Corrector
    if(UseRecoil) {
        recoilCorr = new RecoilCorrector(RecoilMCFile);
        recoilCorr->addDataFile(RecoilDataFile);
        recoilCorr->addMCFile(RecoilMCFile);
    }
    
    
    std::cout << " --- JetAnalyzer initialization ---" << std::endl;
    std::cout << "  jet Id            :\t" << JetId << std::endl;
    std::cout << "  jet pT [1, 2]     :\t" << Jet1Pt << "\t" << Jet2Pt << std::endl;
    std::cout << "  jet eta           :\t" << JetEta << std::endl;
    std::cout << "  b-tagging algo    :\t" << BTag << std::endl;
    std::cout << "  b-tag cut [1, 2]  :\t" << Jet1BTag << "\t" << Jet2BTag << std::endl;
    std::cout << "  apply recoil corr :\t" << (UseRecoil ? "YES" : "NO") << std::endl;
    std::cout << "  recoil file MC    :\t" << RecoilMCFile << std::endl;
    std::cout << "  recoil file Data  :\t" << RecoilDataFile << std::endl;
    std::cout << std::endl;
}

JetAnalyzer::~JetAnalyzer() {
//    JESFile->Close();
    
    if(UseRecoil) delete recoilCorr;
}





std::vector<pat::Jet> JetAnalyzer::FillJetVector(const edm::Event& iEvent) {
    bool isMC(!iEvent.isRealData());
    int BTagTh(Jet1BTag);
    float PtTh(Jet1Pt), EtaTh(JetEta);
    std::vector<pat::Jet> Vect;
    // Declare and open collection
    edm::Handle<std::vector<pat::Jet> > PFJetsCollection;
    iEvent.getByToken(JetToken, PFJetsCollection);
    // Jet corrector
    edm::Handle<reco::JetCorrector> Corrector;
    if(RecalibrateJets || RecalibrateMass) iEvent.getByToken(CorToken, Corrector);
    
    // Open QG value maps
    edm::Handle<edm::ValueMap<float>> QGHandle;
    if(AddQG) iEvent.getByToken(QGToken, QGHandle);
    
    // Loop on Jet collection
    for(std::vector<pat::Jet>::const_iterator it=PFJetsCollection->begin(); it!=PFJetsCollection->end(); ++it) {
        if(Vect.size()>0) {
            PtTh=Jet2Pt;
            BTagTh=Jet2BTag;
        }
        pat::Jet jet=*it;
        int idx=it-PFJetsCollection->begin();
        jet.addUserInt("Index", idx);
        pat::JetRef jetRef(PFJetsCollection, idx);
        // Jet Energy Smearing
        if(isMC) {
            const reco::GenJet* genJet=jet.genJet();
            if(genJet) {
                float smearFactor=GetResolutionRatio(jet.eta());
                reco::Candidate::LorentzVector smearedP4;
                smearedP4=jet.p4()-genJet->p4();
                smearedP4*=smearFactor; // +- 3*smearFactorErr;
                smearedP4+=genJet->p4();
                jet.setP4(smearedP4);
            }
        }
        // Pt and eta cut
        if(jet.pt()<PtTh || fabs(jet.eta())>EtaTh) continue;
        // Quality cut
        if(JetId==1 && !isLooseJet(jet)) continue;
        if(JetId==2 && !isTightJet(jet)) continue;
        if(JetId==3 && !isTightLepVetoJet(jet)) continue;
        // b-tagging
        if(BTagTh==1 && jet.bDiscriminator(BTag)<BTagTh) continue;
        // Fill jet scale uncertainty
        jet.addUserInt("isLoose", isLooseJet(jet) ? 1 : 0);
        jet.addUserInt("isTight", isTightJet(jet) ? 1 : 0);
        jet.addUserInt("isTightLepVeto", isTightLepVetoJet(jet) ? 1 : 0);
        jet.addUserFloat("JESUncertainty", jet.pt()*GetScaleUncertainty(jet));
        jet.addUserFloat("ReshapedDiscriminator", ReshapeBtagDiscriminator(jet)[0]);
        jet.addUserFloat("ReshapedDiscriminatorUp", ReshapeBtagDiscriminator(jet)[1]);
        jet.addUserFloat("ReshapedDiscriminatorDown", ReshapeBtagDiscriminator(jet)[2]);

        // PUPPI soft drop mass for AK8 jets
        if(jet.hasSubjets("SoftDropPuppi")) {
            TLorentzVector puppiSoftdrop, puppiSoftdropSubjet;
            auto const & sdSubjetsPuppi = jet.subjets("SoftDropPuppi");
            for (auto const & it : sdSubjetsPuppi) {
                puppiSoftdropSubjet.SetPtEtaPhiM(it->pt(), it->eta(), it->phi(), it->mass());
                puppiSoftdrop += puppiSoftdropSubjet;
            }
            jet.addUserFloat("ak8PFJetsCHSSoftDropPuppiMass", puppiSoftdrop.M());
        }
        
        // CSV reshaping for soft drop subjets
        if(jet.hasSubjets("SoftDrop")) {
            auto const & sdSubjets = jet.subjets("SoftDrop");
            short nsj = 1;
            for (auto const & it : sdSubjets) {
                pat::Jet subjet = it;
                jet.addUserFloat(Form("ReshapedDiscriminator%d",nsj), ReshapeBtagDiscriminator(subjet)[0]);
                jet.addUserFloat(Form("ReshapedDiscriminatorUp%d",nsj), ReshapeBtagDiscriminator(subjet)[1]);
                jet.addUserFloat(Form("ReshapedDiscriminatorDown%d",nsj), ReshapeBtagDiscriminator(subjet)[2]);
                ++nsj;
            }
        }
                
        // Mass correction
        if(RecalibrateMass) {
            float jec = Corrector->correction(jet);
            if(jet.hasUserFloat("ak8PFJetsCHSPrunedMass")) jet.addUserFloat("ak8PFJetsCHSPrunedMassCorr", jet.userFloat("ak8PFJetsCHSPrunedMass") * jec);
            if(jet.hasUserFloat("ak8PFJetsCHSSoftDropMass")) jet.addUserFloat("ak8PFJetsCHSSoftDropMassCorr", jet.userFloat("ak8PFJetsCHSSoftDropMass") * jec);
            if(jet.hasUserFloat("ak8PFJetsCHSSoftDropPuppiMass")) jet.addUserFloat("ak8PFJetsCHSSoftDropPuppiMassCorr", jet.userFloat("ak8PFJetsCHSSoftDropPuppiMass") * jec);
        }
        //QG tagger for AK4 jets
        if(AddQG && jet.nSubjetCollections()<=0) {
            jet.addUserFloat("QGLikelihood", (*QGHandle)[jetRef]);
        }
        
        Vect.push_back(jet); // Fill vector
    }
    return Vect;
}


void JetAnalyzer::CleanJetsFromMuons(std::vector<pat::Jet>& Jets, std::vector<pat::Muon>& Muons, float angle) {
    for(unsigned int m = 0; m < Muons.size(); m++) {
        for(unsigned int j = 0; j < Jets.size(); ) {
            if(deltaR(Jets[j], Muons[m]) < angle) Jets.erase(Jets.begin() + j);
            else j++;
        }
    }
}

void JetAnalyzer::CleanJetsFromElectrons(std::vector<pat::Jet>& Jets, std::vector<pat::Electron>& Electrons, float angle) {
    for(unsigned int e = 0; e < Electrons.size(); e++) {
        for(unsigned int j = 0; j < Jets.size(); ) {
            if(deltaR(Jets[j], Electrons[e]) < angle) Jets.erase(Jets.begin() + j);
            else j++;
        }
    }
}

int JetAnalyzer::GetNBJets(std::vector<pat::Jet>& Jets) {
    int n(0);
    for(unsigned int i = 0; i < Jets.size(); i++) if(abs(Jets[i].hadronFlavour()) == 5) n++;
    return n;
}

pat::MET JetAnalyzer::FillMetVector(const edm::Event& iEvent) {
    
    edm::Handle<std::vector<pat::MET> > MetCollection;
    iEvent.getByToken(MetToken, MetCollection);
    pat::MET MEt = MetCollection->front();
    MEt.addUserFloat("ptRaw", MEt.uncorPt());
    MEt.addUserFloat("phiRaw", MEt.uncorPhi());
    MEt.addUserFloat("ptType1", MEt.pt());
    MEt.addUserFloat("phiType1", MEt.phi());
    return MEt;
}

void JetAnalyzer::ApplyRecoilCorrections(pat::MET& MET, const reco::Candidate::LorentzVector* GenV, const reco::Candidate::LorentzVector* RecoV, int nJets) {
    double MetPt(MET.pt()), MetPhi(MET.phi()), MetPtScaleUp(MET.pt()), MetPhiScaleUp(MET.phi()), MetPtScaleDown(MET.pt()), MetPhiScaleDown(MET.phi()), MetPtResUp(MET.pt()), MetPhiResUp(MET.phi()), MetPtResDown(MET.pt()), MetPhiResDown(MET.phi());
    double GenPt(0.), GenPhi(0.), LepPt(0.), LepPhi(0.), LepPx(0.), LepPy(0.);
    double RecoilX(0.), RecoilY(0.), Upara(0.), Uperp(0.);
    
    if(GenV) {
        GenPt = GenV->pt();
        GenPhi = GenV->phi();
    }
    else {
        throw cms::Exception("JetAnalyzer", "GenV boson is null. No Recoil Correction can be derived");
        return;
    }
    
    if(RecoV) {
        LepPt = RecoV->pt();
        LepPhi = RecoV->phi();
        LepPx = RecoV->px();
        LepPy = RecoV->py();
        RecoilX = - MET.px() - LepPx;
        RecoilY = - MET.py() - LepPy;
        Upara = (RecoilX*LepPx + RecoilY*LepPy) / LepPt;
        Uperp = (RecoilX*LepPy - RecoilY*LepPx) / LepPt;
    }
    
    // Apply Recoil Corrections
    if(UseRecoil) {
        recoilCorr->CorrectType2(MetPt,          MetPhi,          GenPt, GenPhi, LepPt, LepPhi, Upara, Uperp,  0,  0, nJets);
        recoilCorr->CorrectType2(MetPtScaleUp,   MetPhiScaleUp,   GenPt, GenPhi, LepPt, LepPhi, Upara, Uperp,  3,  0, nJets);
        recoilCorr->CorrectType2(MetPtScaleDown, MetPhiScaleDown, GenPt, GenPhi, LepPt, LepPhi, Upara, Uperp, -3,  0, nJets);
        recoilCorr->CorrectType2(MetPtResUp,     MetPhiResUp,     GenPt, GenPhi, LepPt, LepPhi, Upara, Uperp,  0,  3, nJets);
        recoilCorr->CorrectType2(MetPtResDown,   MetPhiResDown,   GenPt, GenPhi, LepPt, LepPhi, Upara, Uperp,  0, -3, nJets);
    }
    
    // Set userFloats for systematics
    MET.addUserFloat("ptScaleUp", MetPtScaleUp);
    MET.addUserFloat("ptScaleDown", MetPtScaleDown);
    MET.addUserFloat("ptResUp", MetPtResUp);
    MET.addUserFloat("ptResDown", MetPtResDown);
    
    // Set new P4
    MET.setP4(reco::Candidate::PolarLorentzVector(MetPt, MET.eta(), MetPhi, MET.mass()));
}

float JetAnalyzer::GetScaleUncertainty(pat::Jet& jet) {
    if(!isJESFile) return 1.;
    float pt(jet.pt()), eta(fabs(jet.eta()));
    return hist->GetBinContent(hist->FindBin(pt, eta));
}

// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
float JetAnalyzer::GetResolutionRatio(float eta) {
    eta=fabs(eta);
    if(eta>=0.0 && eta<0.5) return 1.052; // +-0.012+0.062-0.061
    if(eta>=0.5 && eta<1.1) return 1.057; // +-0.012+0.056-0.055
    if(eta>=1.1 && eta<1.7) return 1.096; // +-0.017+0.063-0.062
    if(eta>=1.7 && eta<2.3) return 1.134; // +-0.035+0.087-0.085
    if(eta>=2.3 && eta<5.0) return 1.288; // +-0.127+0.155-0.153
    return -1.;
}
float JetAnalyzer::GetResolutionErrorUp(float eta) {
    eta=fabs(eta);
    if(eta>=0.0 && eta<0.5) return 1.115;
    if(eta>=0.5 && eta<1.1) return 1.114;
    if(eta>=1.1 && eta<1.7) return 1.161;
    if(eta>=1.7 && eta<2.3) return 1.228;
    if(eta>=2.3 && eta<5.0) return 1.488;
    return -1.;
}
float JetAnalyzer::GetResolutionErrorDown(float eta) {
    eta=fabs(eta);
    if(eta>=0.0 && eta<0.5) return 0.990;
    if(eta>=0.5 && eta<1.1) return 1.001;
    if(eta>=1.1 && eta<1.7) return 1.032;
    if(eta>=1.7 && eta<2.3) return 1.042;
    if(eta>=2.3 && eta<5.0) return 1.089;
    return -1.;
}

// PFJet Quality ID 2015-2016: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
bool JetAnalyzer::isLooseJet(pat::Jet& jet) {
    if(fabs(jet.eta())<=3.){
        if(jet.neutralHadronEnergyFraction()>=0.99) return false;
        if(jet.neutralEmEnergyFraction()>=0.99) return false;
        if(jet.numberOfDaughters()<=1) return false;
        if(fabs(jet.eta())<=2.4) {
            if(jet.chargedHadronEnergyFraction()<=0.) return false;
            if(jet.chargedEmEnergyFraction()>=0.99) return false;
            if(jet.chargedMultiplicity()<=0) return false;
        }
    }
    else{
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if(jet.neutralMultiplicity()<=10) return false;
    }
    return true;
}

bool JetAnalyzer::isTightJet(pat::Jet& jet) {
    if(fabs(jet.eta())<=3.){
        if(jet.neutralHadronEnergyFraction()>=0.90) return false;
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if(jet.numberOfDaughters()<=1) return false;
        if(fabs(jet.eta())<=2.4) {
            if(jet.chargedHadronEnergyFraction()<=0.) return false;
            if(jet.chargedEmEnergyFraction()>=0.99) return false;
            if(jet.chargedMultiplicity()<=0) return false;
        }
    }
    else{
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if(jet.neutralMultiplicity()<=10) return false;
    }
    return true;
}

bool JetAnalyzer::isTightLepVetoJet(pat::Jet& jet) {
    if(fabs(jet.eta())<=3.){
        if(jet.neutralHadronEnergyFraction()>=0.90) return false;
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if(jet.numberOfDaughters()<=1) return false;
        if(jet.muonEnergyFraction()>=0.8) return false;
        if(fabs(jet.eta())<=2.4) {
            if(jet.chargedHadronEnergyFraction()<=0.) return false;
            if(jet.chargedEmEnergyFraction()>=0.90) return false;
            if(jet.chargedMultiplicity()<=0) return false;
        }
    }
    else{
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if(jet.neutralMultiplicity()<=10) return false;
    }
    return true;
}

std::vector<float> JetAnalyzer::ReshapeBtagDiscriminator(pat::Jet& jet) {
    float pt(jet.pt()), eta(fabs(jet.eta())), discr(jet.bDiscriminator(BTag));
    int hadronFlavour_ = std::abs(jet.hadronFlavour());
    std::vector<float> reshapedDiscr = {discr,discr,discr};
    
    if(UseReshape) {
        BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;
        if (hadronFlavour_ == 5) jf = BTagEntry::FLAV_B;
        else if (hadronFlavour_ == 4) jf = BTagEntry::FLAV_C;
        else if (hadronFlavour_ == 0) jf = BTagEntry::FLAV_UDSG;

        reshapedDiscr[0] = discr*reader->eval(jf, eta, pt); 
        reshapedDiscr[1] = discr*reader_up_jes->eval(jf, eta, pt); 
        reshapedDiscr[2] = discr*reader_down_jes->eval(jf, eta, pt); 
    
//     float reshapedDiscr_up_jes        = discr*reader_up_jes->eval(jf, eta, pt); 
//     float reshapedDiscr_up_lf         = discr*reader_up_lf->eval(jf, eta, pt); 
//     float reshapedDiscr_up_hfstats1   = discr*reader_up_hfstats1->eval(jf, eta, pt); 
//     float reshapedDiscr_up_hfstats2   = discr*reader_up_hfstats2->eval(jf, eta, pt); 
//     float reshapedDiscr_up_cferr1     = discr*reader_up_cferr1->eval(jf, eta, pt); 
//     float reshapedDiscr_up_cferr2     = discr*reader_up_cferr2->eval(jf, eta, pt); 
//     float reshapedDiscr_down_jes        = discr*reader_down_jes->eval(jf, eta, pt); 
//     float reshapedDiscr_down_lf         = discr*reader_down_lf->eval(jf, eta, pt); 
//     float reshapedDiscr_down_hfstats1   = discr*reader_down_hfstats1->eval(jf, eta, pt); 
//     float reshapedDiscr_down_hfstats2   = discr*reader_down_hfstats2->eval(jf, eta, pt); 
//     float reshapedDiscr_down_cferr1     = discr*reader_down_cferr1->eval(jf, eta, pt); 
//     float reshapedDiscr_down_cferr2     = discr*reader_down_cferr2->eval(jf, eta, pt); 
    
//     std::cout << Form("pt, eta, b-tag, flav, reshapedDiscr      : %f, %f, %f, %d, %f\n", pt, eta, discr, jf, reshapedDiscr);
//     std::cout << Form(" jes     : %f / %f\n", reshapedDiscr_up_jes, reshapedDiscr_down_jes);
//     std::cout << Form(" lf      : %f / %f\n", reshapedDiscr_up_lf, reshapedDiscr_down_lf);
//     std::cout << Form(" hfstats1: %f / %f\n", reshapedDiscr_up_hfstats1, reshapedDiscr_down_hfstats1);
//     std::cout << Form(" hfstats2: %f / %f\n", reshapedDiscr_up_hfstats2, reshapedDiscr_down_hfstats2);
//     std::cout << Form(" cferr1  : %f / %f\n", reshapedDiscr_up_cferr1, reshapedDiscr_down_cferr1);
//     std::cout << Form(" cferr2  : %f / %f\n", reshapedDiscr_up_cferr2, reshapedDiscr_down_cferr2);
    }
    return reshapedDiscr;
}


/*
bool JetAnalyzer::isMediumJet(pat::Jet& jet) {
    if(jet.neutralHadronEnergyFraction()>0.95) return false;
    if(jet.neutralEmEnergyFraction()>0.95) return false;
    if(jet.numberOfDaughters()<=1) return false;
    if(fabs(jet.eta())<2.4) {
      if(jet.chargedHadronEnergyFraction()<=0.) return false;
      if(jet.chargedEmEnergyFraction()>0.99) return false;
      if(jet.chargedMultiplicity()<=0) return false;
    }
    return true;
}
bool JetAnalyzer::isTightJet(pat::Jet& jet) {
    if(jet.neutralHadronEnergyFraction()>0.90) return false;
    if(jet.neutralEmEnergyFraction()>0.90) return false;
    if(jet.numberOfDaughters()<=1) return false;
    if(fabs(jet.eta())<2.4) {
      if(jet.chargedHadronEnergyFraction()<=0.) return false;
      if(jet.chargedEmEnergyFraction()>0.99) return false;
      if(jet.chargedMultiplicity()<=0) return false;
    }
    return true;
}

*/


