#include "ElectronAnalyzer.h"



ElectronAnalyzer::ElectronAnalyzer(const edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    ElectronToken(CColl.consumes<std::vector<pat::Electron> >(PSet.getParameter<edm::InputTag>("electrons"))),
    VertexToken(CColl.consumes<reco::VertexCollection>(PSet.getParameter<edm::InputTag>("vertices"))),
    EleVetoIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleVetoIdMap"))),
    EleLooseIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleLooseIdMap"))),
    EleMediumIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleMediumIdMap"))),
    EleTightIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleTightIdMap"))),
    EleHEEPIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleHEEPIdMap"))),
    EleMVANonTrigMediumIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleMVANonTrigMediumIdMap"))),
    EleMVANonTrigTightIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleMVANonTrigTightIdMap"))),
    EleMVATrigMediumIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleMVATrigMediumIdMap"))),
    EleMVATrigTightIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleMVATrigTightIdMap"))),
    Electron1Id(PSet.getParameter<int>("electron1id")),
    Electron2Id(PSet.getParameter<int>("electron2id")),
    Electron1Iso(PSet.getParameter<int>("electron1iso")),
    Electron2Iso(PSet.getParameter<int>("electron2iso")),
    Electron1Pt(PSet.getParameter<double>("electron1pt")),
    Electron2Pt(PSet.getParameter<double>("electron2pt"))
{
    isEleTriggerFile=isEleEfficiencyFile=false;
    
    // AN-13-022
    // Electron trigger
    EleTriggerFile=new TFile("data/DETrigger.root", "READ");
    if(!EleTriggerFile->IsZombie()) {
        EleTriggerDATAHighLeg=(TH2F*)EleTriggerFile->Get("test/DATA_Ele17Leg");
        EleTriggerDATALowLeg=(TH2F*)EleTriggerFile->Get("test/DATA_Ele8Leg");
        EleTriggerMCHighLeg=(TH2F*)EleTriggerFile->Get("test/MC_Ele17Leg");
        EleTriggerMCLowLeg=(TH2F*)EleTriggerFile->Get("test/MC_Ele8Leg");
        isEleTriggerFile=true;
    }
    else std::cout << " - ElectronAnalyzer Warning: No EleTrigger Weight File" << std::endl;
    
    // Electron reco+id+iso
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentification#Triggering_MVA
    EleEfficiencyFile=new TFile("data/electrons_scale_factors.root", "READ");
    if(!EleEfficiencyFile->IsZombie()) {
        ElectronId=(TH2F*)EleEfficiencyFile->Get("electronsDATAMCratio_FO_ID");
        ElectronIso=(TH2F*)EleEfficiencyFile->Get("electronsDATAMCratio_FO_ISO");
        isEleEfficiencyFile=true;
    }
    else std::cout << " - ElectronAnalyzer Warning: No EleEfficiency Weight File" << std::endl;
    
    std::cout << " - ElectronAnalyzer initialized:" << std::endl;
    std::cout << "Id  :\t" << Electron1Id << "\t" << Electron2Id << std::endl;
    std::cout << "Iso :\t" << Electron1Iso << "\t" << Electron2Iso << std::endl;
    std::cout << "pT  :\t" << Electron1Pt << "\t" << Electron2Pt << std::endl;
}

ElectronAnalyzer::~ElectronAnalyzer() {
    EleTriggerFile->Close();
    EleEfficiencyFile->Close();
}





std::vector<pat::Electron> ElectronAnalyzer::FillElectronVector(const edm::Event& iEvent) {
    bool isMC(!iEvent.isRealData());
    int IdTh(Electron1Id), IsoTh(Electron1Iso);
    float PtTh(Electron1Pt);
    std::vector<pat::Electron> Vect;
    // Declare and open collection
    edm::Handle<std::vector<pat::Electron> > EleCollection;
    iEvent.getByToken(ElectronToken, EleCollection);
    
    //edm::Handle<std::vector<pat::Conversion> > EleConv;
    //iEvent.getByToken(edm::InputTag("patConversions"), EleConv);
    
    edm::Handle<reco::VertexCollection> PVCollection;
    iEvent.getByToken(VertexToken, PVCollection);
    const reco::Vertex* vertex=&PVCollection->front();
    
    //value map for ID 2015-2016
    edm::Handle<edm::ValueMap<bool> > VetoIdDecisions;
    edm::Handle<edm::ValueMap<bool> > LooseIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MediumIdDecisions;
    edm::Handle<edm::ValueMap<bool> > TightIdDecisions;
    edm::Handle<edm::ValueMap<bool> > HEEPIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MVANonTrigMediumIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MVANonTrigTightIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MVATrigMediumIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MVATrigTightIdDecisions;
    iEvent.getByToken(EleVetoIdMapToken, VetoIdDecisions);
    iEvent.getByToken(EleLooseIdMapToken, LooseIdDecisions);
    iEvent.getByToken(EleMediumIdMapToken, MediumIdDecisions);
    iEvent.getByToken(EleTightIdMapToken, TightIdDecisions);
    iEvent.getByToken(EleHEEPIdMapToken, HEEPIdDecisions);
    iEvent.getByToken(EleMVANonTrigMediumIdMapToken, MVANonTrigMediumIdDecisions);
    iEvent.getByToken(EleMVANonTrigTightIdMapToken, MVANonTrigTightIdDecisions);
    iEvent.getByToken(EleMVATrigMediumIdMapToken, MVATrigMediumIdDecisions);
    iEvent.getByToken(EleMVATrigTightIdMapToken, MVATrigTightIdDecisions);
    unsigned int elIdx = 0;

    // Loop on Electron collection
    for(std::vector<pat::Electron>::const_iterator it=EleCollection->begin(); it!=EleCollection->end(); ++it) {
        if(Vect.size()>0) {
            IdTh=Electron2Id;
            IsoTh=Electron2Iso;
            PtTh=Electron2Pt;
        }
        pat::Electron el=*it;
	pat::ElectronRef elRef(EleCollection, elIdx);
        // Pt and eta
        if(el.pt()<PtTh || fabs(el.eta())>2.5) continue;
        // PF (?) Isolation R=0.4 https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPFBasedIsolation#for_PAT_electron_users_using_sta
        float pfIso04 = ( el.chargedHadronIso() + std::max(el.neutralHadronIso() + el.photonIso() - 0.5*el.puChargedHadronIso(), 0.) ) / el.pt();
        float pfIso03 = ( el.pfIsolationVariables().sumChargedHadronPt + std::max(el.pfIsolationVariables().sumNeutralHadronEt + el.pfIsolationVariables().sumPhotonEt - 0.5*el.pfIsolationVariables().sumPUPt, 0.) ) / el.pt();
        if(IsoTh==1 && pfIso03>0.15) continue;

        //Electron CutBased and HEEP ID 2015-2016, https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
        bool isPassVeto = (*VetoIdDecisions)[elRef];
        bool isPassLoose = (*LooseIdDecisions)[elRef];
        bool isPassMedium = (*MediumIdDecisions)[elRef];
        bool isPassTight = (*TightIdDecisions)[elRef];
        bool isPassHEEP = (*HEEPIdDecisions)[elRef];
        bool isPassMVANonTrigMedium = (*MVANonTrigMediumIdDecisions)[elRef];
        bool isPassMVANonTrigTight = (*MVANonTrigTightIdDecisions)[elRef];
        bool isPassMVATrigMedium = (*MVATrigMediumIdDecisions)[elRef];
        bool isPassMVATrigTight = (*MVATrigTightIdDecisions)[elRef];

        if(IdTh==0 && !isPassVeto) continue;
        if(IdTh==1 && !isPassLoose) continue;
        if(IdTh==2 && !isPassMedium) continue;
        if(IdTh==3 && !isPassTight) continue;
        if(IdTh==4 && !isPassHEEP) continue;
        if(IdTh==5 && !isPassMVANonTrigMedium) continue;
        if(IdTh==6 && !isPassMVANonTrigTight) continue;
        if(IdTh==7 && !isPassMVATrigMedium) continue;
        if(IdTh==8 && !isPassMVATrigTight) continue;

        ++elIdx;
        
        // Fill userFloat
        el.addUserFloat("pfIso03", pfIso03);
        el.addUserFloat("pfIso04", pfIso04);
        el.addUserFloat("trkIso", el.pfIsolationVariables().sumChargedHadronPt);
        el.addUserFloat("dxy", el.gsfTrack()->dxy( vertex->position() ));
        el.addUserFloat("dz", el.gsfTrack()->dz( vertex->position() ));
        el.addUserInt("isVeto", isPassVeto ? 1 : 0);
        el.addUserInt("isLoose", isPassLoose ? 1 : 0);
        el.addUserInt("isMedium", isPassMedium ? 1 : 0);
        el.addUserInt("isTight", isPassTight ? 1 : 0);
        el.addUserInt("isHEEP", isPassHEEP ? 1 : 0);
        el.addUserInt("isMVANonTrigMedium", isPassMVANonTrigMedium ? 1 : 0);
        el.addUserInt("isMVANonTrigTight", isPassMVANonTrigTight ? 1 : 0);
        el.addUserInt("isMVATrigMedium", isPassMVATrigMedium ? 1 : 0);
        el.addUserInt("isMVATrigTight", isPassMVATrigTight ? 1 : 0);
        // Fill vector
        Vect.push_back(el);
    }
    return Vect;
}

/*

//// Obsolete Electron ID for Run1

// Electron Cut-based Quality ID: see https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification#Electron_ID_Working_Points
bool ElectronAnalyzer::IsElectronCut(pat::Electron& el, const reco::Vertex* vertex, int q) {
    // Electron Quality: 0 - Veto, 1 - Loose, 2 - Medium, 3 - Tight
    //                                VBTF90     VBTF 80    VBTF 70
    if(q<0 || q>3) return false;
    //     ECAL BARREL                                   ECAL ENDCAP
    float sieieEB[4]    ={0.01,  0.01,  0.01,  0.01 },  sieieEE[4]    ={0.03,  0.03,  0.03,  0.03 };
    float dEtaInEB[4]   ={0.007, 0.007, 0.004, 0.004},  dEtaInEE[4]   ={0.010, 0.009, 0.007, 0.005};
    float dPhiInEB[4]   ={0.8,   0.15,  0.06,  0.03 },  dPhiInEE[4]   ={0.7,   0.10,  0.03,  0.02 };
    float hOeEB[4]      ={0.15,  0.12,  0.12,  0.12 },  hOeEE[4]      ={1.e10, 0.10,  0.10,  0.10 };
    float d0vtxEB[4]    ={0.04,  0.02,  0.02,  0.02 },  d0vtxEE[4]    ={0.04,  0.02,  0.02,  0.02 };
    float dZvtxEB[4]    ={0.2,   0.2,   0.1,   0.1  },  dZvtxEE[4]    ={0.2,   0.2,   0.1,   0.1  };
    float abs1oE1opEB[4]={0.,    0.05,  0.05,  0.05 },  abs1oE1opEE[4]={0.,    0.05,  0.05,  0.05 };
    //float pfIsoOPtEB[4] ={0.15,  0.15,  0.15,  0.10 },  pfIsoOPtEE[4] ={0.15,  0.15,  0.15,  0.10 };
    //float vtxFitEB[4]   ={1.e-6, 1.e-6, 1.e-6, 1.e-6},  vtxFitEE[4]   ={1.e-6, 1.e-6, 1.e-6, 1.e-6};
    float mHitsEB[4]    ={1.,    1.,    1.,    0.   },  mHitsEE[4]    ={1.,    1.,    1.,    0.   };
    
    bool isEB=el.isEB();
    bool isEE=el.isEE();
    
    float sieie=el.scSigmaIEtaIEta();
    float dEtaIn=fabs(el.deltaEtaSuperClusterTrackAtVtx());
    float dPhiIn=fabs(el.deltaPhiSuperClusterTrackAtVtx());
    float hOe=el.hadronicOverEm();
    float d0vtx=fabs( el.gsfTrack()->d0() - vertex->x()*sin(el.gsfTrack()->phi()) + vertex->y()*cos(el.gsfTrack()->phi()) );
    float dZvtx=fabs( (el.vz()-vertex->z()) -( (el.vx()-vertex->x())*el.px()+(el.vy()-vertex->y())*el.py() )/el.pt()*el.pz()/el.pt() );
    float E = el.superCluster()->energy();
    float p = 1./(el.eSuperClusterOverP()/E);
    float abs1oE1op=fabs(1./E - 1./p);
    float mHits=el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    
    if(isEB) if(sieie>sieieEB[q] || dEtaIn>dEtaInEB[q] || dPhiIn>dPhiInEB[q] || hOe>hOeEB[q] || d0vtx>d0vtxEB[q] || dZvtx>dZvtxEB[q] || abs1oE1op>abs1oE1opEB[q] || mHits>mHitsEB[q]) return false;
    if(isEE) if(sieie>sieieEE[q] || dEtaIn>dEtaInEE[q] || dPhiIn>dPhiInEE[q] || hOe>hOeEE[q] || d0vtx>d0vtxEE[q] || dZvtx>dZvtxEE[q] || abs1oE1op>abs1oE1opEE[q] || mHits>mHitsEE[q]) return false;
    return true;
}

// https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification#Glossary_of_Variables
bool ElectronAnalyzer::IsLooseElectron(pat::Electron& el, const reco::Vertex* vertex) {
    if(el.isEB()) {
        if(  
               fabs(el.deltaEtaSuperClusterTrackAtVtx())<0.007
            && fabs(el.deltaPhiSuperClusterTrackAtVtx())<0.15
            && el.scSigmaIEtaIEta()<0.01
            && el.hadronicOverEm()<0.12
            && el.gsfTrack()->dxy(vertex->position())<0.02
            && el.gsfTrack()->dz(vertex->position())<0.2
            && fabs(1./el.ecalEnergy() - 1./el.trackMomentumAtVtx().R())<0.05
      //      && el.reco::GsfElectron::fbrem()<10.e-6
            && el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=1
        ) return true;
    }
    else if(el.isEE()) {
        if(
               fabs(el.deltaEtaSuperClusterTrackAtVtx())<0.009
            && fabs(el.deltaPhiSuperClusterTrackAtVtx())<0.10
            && el.scSigmaIEtaIEta()<0.03
            && el.hadronicOverEm()<0.10
            && el.gsfTrack()->dxy(vertex->position())<0.02
            && el.gsfTrack()->dz(vertex->position())<0.2
            && fabs(1./el.ecalEnergy() - 1./el.trackMomentumAtVtx().R())<0.05
      //      && el.reco::GsfElectron::fbrem()<10.e-6
            && el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=1
        ) return true;
    }
    return false;
}



bool ElectronAnalyzer::IsLooseMVAElectron(pat::Electron& el) {
    //if(el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)>0) return false;
    float mva=el.electronID("mvaTrigV0");
    float eta=fabs(el.eta());
    // AN2013_108_v8 page 31
    if(el.pt()>=10.) {
        if(eta<0.8)               {if(mva<-0.34) return false;}
        if(eta>=0.8 && eta<1.479) {if(mva<-0.65) return false;}
        if(eta>=1.479)            {if(mva<+0.60) return false;}
    }
    else {
        if(eta<0.8)               {if(mva<+0.47) return false;}
        if(eta>=0.8 && eta<1.479) {if(mva<0.004) return false;}
        if(eta>=1.479)            {if(mva<0.295) return false;}
    }
    return true;
}

// https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification#Recommended_Working_Points_With
bool ElectronAnalyzer::IsTightMVAElectron(pat::Electron& el) {
    // Conversion Veto
    if(!el.passConversionVeto()) return false;
    // Missing hits veto
    if(el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)>0) return false;
    // Triggering MVA Cut
    float mva=el.electronID("mvaTrigV0");
    float eta=fabs(el.eta());
    if(el.pt()<20.) {
        if(eta<0.8)               {if(mva<+0.00) return false;}
        if(eta>=0.8 && eta<1.479) {if(mva<+0.10) return false;}
        if(eta>=1.479)            {if(mva<+0.62) return false;}
    }
    else {
        if(eta<0.8)               {if(mva<+0.94) return false;}
        if(eta>=0.8 && eta<1.479) {if(mva<+0.85) return false;}
        if(eta>=1.479)            {if(mva<+0.92) return false;}
    }
    return true;
}

*/

float ElectronAnalyzer::GetDoubleElectronTriggerSF(pat::Electron& el1, pat::Electron& el2) {
    if(!isEleTriggerFile) return 1.;
    float pt1=el1.pt(), pt2=el2.pt();
    float eta1(fabs(el1.eta())), eta2(fabs(el2.eta()));
    if(pt1>=EleTriggerPtMax) pt1=EleTriggerPtMax-1.;
    if(pt2>=EleTriggerPtMax) pt2=EleTriggerPtMax-1.;
    
    float effDATAHighEle1 = EleTriggerDATAHighLeg->GetBinContent(EleTriggerDATAHighLeg->FindBin(pt1, eta1));
    float effDATALowEle2 = EleTriggerDATALowLeg->GetBinContent(EleTriggerDATALowLeg->FindBin(pt2, eta2));
    float effDATAHighEle2 = EleTriggerDATAHighLeg->GetBinContent(EleTriggerDATAHighLeg->FindBin(pt2, eta2));
    float effDATALowEle1 = EleTriggerDATALowLeg->GetBinContent(EleTriggerDATALowLeg->FindBin(pt1, eta1));
    
    float effMCHighEle1 = EleTriggerMCHighLeg->GetBinContent(EleTriggerMCHighLeg->FindBin(pt1, eta1));
    float effMCLowEle2 = EleTriggerMCLowLeg->GetBinContent(EleTriggerMCLowLeg->FindBin(pt2, eta2));
    float effMCHighEle2 = EleTriggerMCHighLeg->GetBinContent(EleTriggerMCHighLeg->FindBin(pt2, eta2));
    float effMCLowEle1 = EleTriggerMCLowLeg->GetBinContent(EleTriggerMCLowLeg->FindBin(pt1, eta1));
    
    float effDATA = effDATAHighEle1*effDATALowEle2 + effDATALowEle1*effDATAHighEle2 - effDATAHighEle1*effDATAHighEle2;
    float effMC = effMCHighEle1*effMCLowEle2 + effMCLowEle1*effMCHighEle2 - effMCHighEle1*effMCHighEle2;
    return effDATA/effMC;
}

float ElectronAnalyzer::GetElectronIdSF(pat::Electron& el) {
    if(!isEleEfficiencyFile) return 1.;
    float pt(el.pt()), eta(fabs(el.eta()));
    if(pt>=EleEfficiencyPtMax) pt=EleEfficiencyPtMax-1.;
    return ElectronId->GetBinContent(ElectronId->FindBin(eta, pt));
}

float ElectronAnalyzer::GetElectronIdSFError(pat::Electron& el) {
    if(!isEleEfficiencyFile) return 1.;
    float pt(el.pt()), eta(fabs(el.eta()));
    if(pt>=EleEfficiencyPtMax) pt=EleEfficiencyPtMax-1.;
    return ElectronId->GetBinError(ElectronId->FindBin(eta, pt));
}

float ElectronAnalyzer::GetElectronIsoSF(pat::Electron& el) {
    if(!isEleEfficiencyFile) return 1.;
    float pt(el.pt()), eta(fabs(el.eta()));
    if(pt>=EleEfficiencyPtMax) pt=EleEfficiencyPtMax-1.;
    return ElectronIso->GetBinContent(ElectronIso->FindBin(eta, pt));
}

float ElectronAnalyzer::GetElectronIsoSFError(pat::Electron& el) {
    if(!isEleEfficiencyFile) return 1.;
    float pt(el.pt()), eta(fabs(el.eta()));
    if(pt>=EleEfficiencyPtMax) pt=EleEfficiencyPtMax-1.;
    return ElectronIso->GetBinError(ElectronIso->FindBin(eta, pt));
}

// https://twiki.cern.ch/twiki/bin/view/Main/EGammaScaleFactors2012#2012_8_TeV_data_53X
float ElectronAnalyzer::GetLooseElectronSF(pat::Electron& el) {
    float pt=el.pt();
    float eta=fabs(el.eta());
    if(pt>10. && pt<=15.) {
        if(eta<=0.8)                     return 0.855;
        else if(eta>0.8   && eta<=1.442) return 0.858;
        else if(eta>1.442 && eta<=1.556) return 1.109;
        else if(eta>1.556 && eta<=2.0  ) return 0.838;
        else if(eta>2.0   && eta<=2.5  ) return 1.034;
    }
    else if(pt>15. && pt<=20.) {
        if(eta<=0.8)                     return 0.962;
        else if(eta>0.8   && eta<=1.442) return 0.962;
        else if(eta>1.442 && eta<=1.556) return 0.903;
        else if(eta>1.556 && eta<=2.0  ) return 0.939;
        else if(eta>2.0   && eta<=2.5  ) return 0.970;
    }
    else if(pt>20. && pt<=30.) {
        if(eta<=0.8)                     return 1.005;
        else if(eta>0.8   && eta<=1.442) return 0.981;
        else if(eta>1.442 && eta<=1.556) return 1.044;
        else if(eta>1.556 && eta<=2.0  ) return 0.980;
        else if(eta>2.0   && eta<=2.5  ) return 1.017;
    }
    else if(pt>30. && pt<=40.) {
        if(eta<=0.8)                     return 1.004;
        else if(eta>0.8   && eta<=1.442) return 0.991;
        else if(eta>1.442 && eta<=1.556) return 0.998;
        else if(eta>1.556 && eta<=2.0  ) return 0.992;
        else if(eta>2.0   && eta<=2.5  ) return 1.019;
    }
    else if(pt>40. && pt<=50.) {
        if(eta<=0.8)                     return 1.008;
        else if(eta>0.8   && eta<=1.442) return 0.994;
        else if(eta>1.442 && eta<=1.556) return 0.989;
        else if(eta>1.556 && eta<=2.0  ) return 1.004;
        else if(eta>2.0   && eta<=2.5  ) return 1.005;
    }
    else if(pt>50.) {
        if(eta<=0.8)                     return 1.008;
        else if(eta>0.8   && eta<=1.442) return 0.999;
        else if(eta>1.442 && eta<=1.556) return 0.994;
        else if(eta>1.556 && eta<=2.0  ) return 1.006;
        else if(eta>2.0   && eta<=2.5  ) return 1.009;
    }
    return 1.;
}


