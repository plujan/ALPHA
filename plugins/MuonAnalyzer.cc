#include "MuonAnalyzer.h"
  
MuonAnalyzer::MuonAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    MuonToken(CColl.consumes<std::vector<pat::Muon> >(PSet.getParameter<edm::InputTag>("muons"))),
    VertexToken(CColl.consumes<reco::VertexCollection>(PSet.getParameter<edm::InputTag>("vertices"))),
    MuonIdFileName(PSet.getParameter<std::string>("muonIdFileName")),
    MuonIsoFileName(PSet.getParameter<std::string>("muonIsoFileName")),
    MuonHighptFileName(PSet.getParameter<std::string>("muonHighptFileName")),
    MuonTriggerFileName(PSet.getParameter<std::string>("muonTriggerFileName")),
    DoubleMuonTriggerFileName(PSet.getParameter<std::string>("doubleMuonTriggerFileName")), //obsolete
    Muon1Id(PSet.getParameter<int>("muon1id")),
    Muon2Id(PSet.getParameter<int>("muon2id")),
    Muon1Iso(PSet.getParameter<int>("muon1iso")),
    Muon2Iso(PSet.getParameter<int>("muon2iso")),
    Muon1Pt(PSet.getParameter<double>("muon1pt")),
    Muon2Pt(PSet.getParameter<double>("muon2pt"))
{
    isMuonTriggerFile = isDoubleMuonTriggerFile = isMuonIdFile = isMuonHighptFile = false;
    
    // Double Muon trigger: obsolete!
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs
    DoubleMuonTriggerFile=new TFile(DoubleMuonTriggerFileName.c_str(), "READ");
    if(!DoubleMuonTriggerFile->IsZombie()) {
        MuonTriggerLt20=(TH2F*)DoubleMuonTriggerFile->Get("DATA_over_MC_Mu17Mu8_Tight_Mu1_10To20_&_Mu2_20ToInfty_with_SYST_uncrt");
        MuonTriggerGt20=(TH2F*)DoubleMuonTriggerFile->Get("DATA_over_MC_Mu17Mu8_Tight_Mu1_20ToInfty_&_Mu2_20ToInfty_with_SYST_uncrt");
        for(int i=1; i<=MuonTriggerGt20->GetNbinsX(); i++) {
            for(int j=1; j<=MuonTriggerGt20->GetNbinsY(); j++) {
                if(j>i) {
                    if(MuonTriggerGt20->GetBinContent(i, j)>0.) std::cout << " - MuonAnalyzer Warning: Trying to symmetrize diagonal matrix in bin " << i << ", " << j << std::endl;
                    MuonTriggerGt20->SetBinContent(i, j, MuonTriggerGt20->GetBinContent(j, i));
                    MuonTriggerGt20->SetBinError(i, j, MuonTriggerGt20->GetBinError(j, i));
                }
            }
        }
        isDoubleMuonTriggerFile=true;
    }
    else {
        throw cms::Exception("MuonAnalyzer", "No Double Muon Trigger Weight File");
        return;
    }

    //Single Muon Trigger, 2015-2016
    MuonTriggerFile=new TFile(MuonTriggerFileName.c_str(), "READ");
    if(!MuonTriggerFile->IsZombie()) {
        MuonTriggerIsoMu20=(TH2F*)MuonTriggerFile->Get("runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins/pt_abseta_ratio");
        isMuonTriggerFile=true;
    }
    else {
        throw cms::Exception("MuonAnalyzer", "No Muon Trigger Weight File");
        return;
    }
    
    //Muon id and iso, 2015-2016
    MuonIdFile=new TFile(MuonIdFileName.c_str(), "READ");
    if(!MuonIdFile->IsZombie()) {
        MuonIdLoose=(TH2F*)MuonIdFile->Get("NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1/pt_abseta_ratio");
        MuonIdTight=(TH2F*)MuonIdFile->Get("NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/pt_abseta_ratio");
        isMuonIdFile=true;
    }
    else {
        throw cms::Exception("MuonAnalyzer", "No MuonId Weight File");
        return;
    }

    MuonIsoFile=new TFile(MuonIsoFileName.c_str(), "READ");
    if(!MuonIsoFile->IsZombie()) {
        MuonIsoLoose=(TH2F*)MuonIsoFile->Get("NUM_LooseRelIso_DEN_LooseID_PAR_pt_spliteta_bin1/pt_abseta_ratio");
        MuonIsoTight=(TH2F*)MuonIsoFile->Get("NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/pt_abseta_ratio");
        isMuonIsoFile=true;
    }
    else {
        throw cms::Exception("MuonAnalyzer", "No MuonIso Weight File");
        return;
    }

    MuonHighptFile=new TFile(MuonHighptFileName.c_str(), "READ");
    if(!MuonHighptFile->IsZombie()) {
        MuonIdHighpt=(TH2F*)MuonHighptFile->Get("HighPtID_PtEtaBins_Pt53/pTtuneP_abseta_ratio");
        MuonIsoHighpt=(TH2F*)MuonHighptFile->Get("tkRelIsoID_PtEtaBins_Pt53/pTtuneP_abseta_ratio");
        isMuonHighptFile=true;
    }
    else {
        throw cms::Exception("MuonAnalyzer", "No MuonHighpt Weight File");
        return;
    }

}

MuonAnalyzer::~MuonAnalyzer() {
    DoubleMuonTriggerFile->Close();
    MuonTriggerFile->Close();
    MuonIdFile->Close();
    MuonIsoFile->Close();
    MuonHighptFile->Close();
}





std::vector<pat::Muon> MuonAnalyzer::FillMuonVector(const edm::Event& iEvent) {
    bool isMC(!iEvent.isRealData());
    int IdTh(Muon1Id), IsoTh(Muon1Iso);
    float PtTh(Muon1Pt);
    std::vector<pat::Muon> Vect;
    // Declare and open collections
    edm::Handle<std::vector<pat::Muon> > MuonCollection;
    iEvent.getByToken(MuonToken, MuonCollection);
    
    edm::Handle<reco::VertexCollection> PVCollection;
    iEvent.getByToken(VertexToken, PVCollection);
    const reco::Vertex* vertex=&PVCollection->front();
    
    // Loop on Muon collection
    for(std::vector<pat::Muon>::const_iterator it=MuonCollection->begin(); it!=MuonCollection->end(); ++it) {
        if(Vect.size()>0) {
            IdTh=Muon2Id;
            IsoTh=Muon2Iso;
            PtTh=Muon2Pt;
        }
        pat::Muon mu=*it;
        // Pt and eta
        if(mu.pt()<PtTh || fabs(mu.eta())>2.4) continue;
        // Isolation
        float pfIso03 = (mu.pfIsolationR03().sumChargedHadronPt + std::max(mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - 0.5*mu.pfIsolationR03().sumPUPt, 0.) ) / mu.pt();
        float pfIso04 = (mu.pfIsolationR04().sumChargedHadronPt + std::max(mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt, 0.) ) / mu.pt();
        // Muon Isolation working point 2015-2016: see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
        if(IsoTh==1 && pfIso04>0.25) continue;
        if(IsoTh==2 && pfIso04>0.15) continue;
        // Muon Quality ID 2015-2016: see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        if(IdTh==1 && !mu.isLooseMuon()) continue;
        if(IdTh==2 && !mu.isMediumMuon()) continue;
        if(IdTh==3 && !mu.isTightMuon(*vertex)) continue;
        if(IdTh==4 && !mu.isHighPtMuon(*vertex)) continue;
        // Add userFloat
        mu.addUserFloat("pfIso03", pfIso03);
        mu.addUserFloat("pfIso04", pfIso04);
        mu.addUserFloat("dxy", mu.muonBestTrack()->dxy(vertex->position()));
        mu.addUserFloat("dz", mu.muonBestTrack()->dz(vertex->position()));
        mu.addUserInt("isVeto", IsCustomTracker(mu, vertex) ? 1 : 0);
        mu.addUserInt("isLoose", mu.isLooseMuon() ? 1 : 0);
        mu.addUserInt("isMedium", mu.isMediumMuon() ? 1 : 0);
        mu.addUserInt("isTight", mu.isTightMuon(*vertex) ? 1 : 0);
        mu.addUserInt("isHighPt", mu.isHighPtMuon(*vertex) ? 1 : 0);
        // Fill vector
        Vect.push_back(mu);
    }
    return Vect;
}


bool MuonAnalyzer::IsCustomTracker(pat::Muon& mu, const reco::Vertex* vertex) {
    return false;
}


//obsolete
float MuonAnalyzer::GetDoubleMuonTriggerSF(pat::Muon& mu1, pat::Muon& mu2) {
    if(!isDoubleMuonTriggerFile) return 1.;
    float eta1=fabs(mu1.eta());
    float eta2=fabs(mu2.eta());
    // Muon POG enumeration is inverted 1 <-> 2
    if(mu2.pt()<20.) return MuonTriggerLt20->GetBinContent(MuonTriggerLt20->FindBin(eta2, eta1));
    return MuonTriggerGt20->GetBinContent(MuonTriggerGt20->FindBin(eta2, eta1));
}

//obsolete
float MuonAnalyzer::GetDoubleMuonTriggerSFError(pat::Muon& mu1, pat::Muon& mu2) {
    if(!isDoubleMuonTriggerFile) return 1.;
    float eta1=fabs(mu1.eta());
    float eta2=fabs(mu2.eta());
    // Muon POG enumeration is inverted 1 <-> 2
    if(mu2.pt()<20.) return MuonTriggerLt20->GetBinError(MuonTriggerLt20->FindBin(eta2, eta1));
    return MuonTriggerGt20->GetBinError(MuonTriggerGt20->FindBin(eta2, eta1));
}

float MuonAnalyzer::GetMuonTriggerSFIsoMu20(pat::Muon& mu) {
    if(!isMuonTriggerFile) return 1.;
    double pt = std::min( std::max( MuonTriggerIsoMu20->GetXaxis()->GetXmin(), mu.pt() ) , MuonTriggerIsoMu20->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonTriggerIsoMu20->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonTriggerIsoMu20->GetBinContent( MuonTriggerIsoMu20->FindBin(pt, abseta) );
}

float MuonAnalyzer::GetMuonTriggerSFErrorIsoMu20(pat::Muon& mu) {
    if(!isMuonTriggerFile) return 1.;
    double pt = std::min( std::max( MuonTriggerIsoMu20->GetXaxis()->GetXmin(), mu.pt() ) , MuonTriggerIsoMu20->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonTriggerIsoMu20->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonTriggerIsoMu20->GetBinError( MuonTriggerIsoMu20->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonIdSFLoose(pat::Muon& mu) {
    if(!isMuonIdFile) return 1.;
    double pt = std::min( std::max( MuonIdLoose->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdLoose->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIdLoose->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIdLoose->GetBinContent( MuonIdLoose->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonIdSFLooseError(pat::Muon& mu) {
    if(!isMuonIdFile) return 1.;
    double pt = std::min( std::max( MuonIdLoose->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdLoose->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIdLoose->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIdLoose->GetBinError( MuonIdLoose->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonIdSFTight(pat::Muon& mu) {
    if(!isMuonIdFile) return 1.;
    double pt = std::min( std::max( MuonIdTight->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdTight->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIdTight->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIdTight->GetBinContent( MuonIdTight->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonIdSFTightError(pat::Muon& mu) {
    if(!isMuonIdFile) return 1.;
    double pt = std::min( std::max( MuonIdTight->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdTight->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIdTight->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIdTight->GetBinError( MuonIdTight->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonIdSFHighpt(pat::Muon& mu) {
    if(!isMuonHighptFile) return 1.;
    double pt = std::min( std::max( MuonIdHighpt->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdHighpt->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIdHighpt->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIdHighpt->GetBinContent( MuonIdHighpt->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonIdSFHighptError(pat::Muon& mu) {
    if(!isMuonHighptFile) return 1.;
    double pt = std::min( std::max( MuonIdHighpt->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdHighpt->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIdHighpt->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIdHighpt->GetBinError( MuonIdHighpt->FindBin(pt,abseta) );
}


float MuonAnalyzer::GetMuonIsoSFLoose(pat::Muon& mu) {
    if(!isMuonIsoFile) return 1.;
    double pt = std::min( std::max( MuonIsoLoose->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoLoose->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIsoLoose->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIsoLoose->GetBinContent( MuonIsoLoose->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonIsoSFLooseError(pat::Muon& mu) {
    if(!isMuonIsoFile) return 1.;
    double pt = std::min( std::max( MuonIsoLoose->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoLoose->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIsoLoose->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIsoLoose->GetBinError( MuonIsoLoose->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonIsoSFTight(pat::Muon& mu) {
    if(!isMuonIsoFile) return 1.;
    double pt = std::min( std::max( MuonIsoTight->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoTight->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIsoTight->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIsoTight->GetBinContent( MuonIsoTight->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonIsoSFTightError(pat::Muon& mu) {
    if(!isMuonIsoFile) return 1.;
    double pt = std::min( std::max( MuonIsoTight->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoTight->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIsoTight->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIsoTight->GetBinError( MuonIsoTight->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonIsoSFHighpt(pat::Muon& mu) {
    if(!isMuonHighptFile) return 1.;
    double pt = std::min( std::max( MuonIsoHighpt->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoHighpt->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIsoHighpt->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIsoHighpt->GetBinContent( MuonIsoHighpt->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonIsoSFHighptError(pat::Muon& mu) {
    if(!isMuonHighptFile) return 1.;
    double pt = std::min( std::max( MuonIsoHighpt->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoHighpt->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonIsoHighpt->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonIsoHighpt->GetBinError( MuonIsoHighpt->FindBin(pt,abseta) );
}

TH1F* MuonAnalyzer::ConvertTGraph(TGraphAsymmErrors* g) {
    int n=g->GetN();
    float x[n+1];
    for(int i=0; i<n; i++) x[i]=g->GetX()[i]-g->GetEXlow()[i];
    x[n]=g->GetX()[n-1]+g->GetEXhigh()[n-1];
    
    TH1F* h=new TH1F(g->GetName(), g->GetTitle(), n, x); h->Sumw2();
    for(int i=0; i<n; i++) {
      h->SetBinContent(i+1, g->GetY()[i]);
      h->SetBinError(i+1, g->GetEYhigh()[i]);
    }
    return h;
}

