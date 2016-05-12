#include "MuonAnalyzer.h"
  
MuonAnalyzer::MuonAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    MuonToken(CColl.consumes<std::vector<pat::Muon> >(PSet.getParameter<edm::InputTag>("muons"))),
    VertexToken(CColl.consumes<reco::VertexCollection>(PSet.getParameter<edm::InputTag>("vertices"))),
    Muon1Id(PSet.getParameter<int>("muon1id")),
    Muon2Id(PSet.getParameter<int>("muon2id")),
    Muon1Iso(PSet.getParameter<int>("muon1iso")),
    Muon2Iso(PSet.getParameter<int>("muon2iso")),
    Muon1Pt(PSet.getParameter<double>("muon1pt")),
    Muon2Pt(PSet.getParameter<double>("muon2pt"))
{
    isMuonTriggerFile=isMuonIdFile=isMuonIsoFile=false;
    
    // Muon trigger
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs
    MuonTriggerFile=new TFile("data/MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.root", "READ");
    if(!MuonTriggerFile->IsZombie()) {
      MuonTriggerLt20=(TH2F*)MuonTriggerFile->Get("DATA_over_MC_Mu17Mu8_Tight_Mu1_10To20_&_Mu2_20ToInfty_with_SYST_uncrt");
      MuonTriggerGt20=(TH2F*)MuonTriggerFile->Get("DATA_over_MC_Mu17Mu8_Tight_Mu1_20ToInfty_&_Mu2_20ToInfty_with_SYST_uncrt");
      for(int i=1; i<=MuonTriggerGt20->GetNbinsX(); i++) {
        for(int j=1; j<=MuonTriggerGt20->GetNbinsY(); j++) {
          if(j>i) {
            if(MuonTriggerGt20->GetBinContent(i, j)>0.) std::cout << " - MuonAnalyzer Warning: Trying to symmetrize diagonal matrix in bin " << i << ", " << j << std::endl;
            MuonTriggerGt20->SetBinContent(i, j, MuonTriggerGt20->GetBinContent(j, i));
            MuonTriggerGt20->SetBinError(i, j, MuonTriggerGt20->GetBinError(j, i));
          }
        }
      }
      isMuonTriggerFile=true;
    }
    else std::cout << " - MuonAnalyzer Warning: No Muon Trigger Weight File" << std::endl;
    
    // Muon ID and ISO
    // https://indico.cern.ch/getFile.py/access?contribId=1&resId=2&materialId=slides&confId=257630
    MuonIdFile=new TFile("data/MuonEfficiencies_Run2012ReReco_53X.root", "READ");
    if(!MuonIdFile->IsZombie()) {
      //DATA_over_MC_Tight_eta_pt20-500
      MuonIdLt09=ConvertTGraph( (TGraphAsymmErrors*)MuonIdFile->Get("DATA_over_MC_Tight_pt_abseta<0.9") );
      MuonId0912=ConvertTGraph( (TGraphAsymmErrors*)MuonIdFile->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2") );
      MuonId1221=ConvertTGraph( (TGraphAsymmErrors*)MuonIdFile->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1") );
      MuonId2124=ConvertTGraph( (TGraphAsymmErrors*)MuonIdFile->Get("DATA_over_MC_Tight_pt_abseta2.1-2.4") );
      isMuonIdFile=true;
    }
    else std::cout << " - MuonAnalyzer Warning: No MuonId Weight File" << std::endl;
    
    MuonIsoFile=new TFile("data/MuonEfficiencies_ISO_Run_2012ReReco_53X.root", "READ");
    if(!MuonIsoFile->IsZombie()) {
      //DATA_over_MC_combRelIsoPF04dBeta<012_Tight_eta_pt20-500
      MuonIsoLt09=ConvertTGraph( (TGraphAsymmErrors*)MuonIsoFile->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta<0.9") );
      MuonIso0912=ConvertTGraph( (TGraphAsymmErrors*)MuonIsoFile->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta0.9-1.2") );
      MuonIso1221=ConvertTGraph( (TGraphAsymmErrors*)MuonIsoFile->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta1.2-2.1") );
      MuonIso2124=ConvertTGraph( (TGraphAsymmErrors*)MuonIsoFile->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta2.1-2.4") );
      isMuonIsoFile=true;
    }
    else std::cout << " - MuonAnalyzer Warning: No MuonIso Weight File" << std::endl;
    
    MuonPtMax=MuonIdLt09->GetXaxis()->GetXmax();
    MuonPtMin=MuonIdLt09->GetXaxis()->GetXmin();
    
    std::cout << " - MuonAnalyzer initialized:" << std::endl;
    std::cout << "Id  :\t" << Muon1Id << "\t" << Muon2Id << std::endl;
    std::cout << "Iso :\t" << Muon1Iso << "\t" << Muon2Iso << std::endl;
    std::cout << "pT  :\t" << Muon1Pt << "\t" << Muon2Pt << std::endl;
}

MuonAnalyzer::~MuonAnalyzer() {
    MuonTriggerFile->Close();
    MuonIdFile->Close();
    MuonIsoFile->Close();
}





std::vector<pat::Muon> MuonAnalyzer::FillMuonVector(const edm::Event& iEvent) {
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
        if(IsoTh==1 && pfIso04>0.20) continue;
        if(IsoTh==2 && pfIso04>0.12) continue;
        // Muon Quality ID: see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
        if(IdTh==1 && !mu.isLooseMuon()) continue;
        if(IdTh==2 && !mu.isMediumMuon()) continue;
        if(IdTh==3 && !mu.isTightMuon(*vertex)) continue;
        // Add userFloat
        mu.addUserFloat("pfIso03", pfIso03);
        mu.addUserFloat("pfIso04", pfIso04);
        mu.addUserFloat("dxy", mu.muonBestTrack()->dxy(vertex->position()));
        mu.addUserFloat("dz", mu.muonBestTrack()->dz(vertex->position()));
        mu.addUserInt("isVeto", IsCustomTracker(mu, vertex) ? 1 : 0);
        mu.addUserInt("isLoose", mu.isLooseMuon() ? 1 : 0);
        mu.addUserInt("isMedium", mu.isMediumMuon() ? 1 : 0);
        mu.addUserInt("isTight", mu.isTightMuon(*vertex) ? 1 : 0);
        // Fill vector
        Vect.push_back(mu);
    }
    return Vect;
}


bool MuonAnalyzer::IsCustomTracker(pat::Muon& mu, const reco::Vertex* vertex) {
    return false;
}


float MuonAnalyzer::GetDoubleMuonTriggerSF(pat::Muon& mu1, pat::Muon& mu2) {
    if(!isMuonTriggerFile) return 1.;
    float eta1=fabs(mu1.eta());
    float eta2=fabs(mu2.eta());
    // Muon POG enumeration is inverted 1 <-> 2
    if(mu2.pt()<20.) return MuonTriggerLt20->GetBinContent(MuonTriggerLt20->FindBin(eta2, eta1));
    return MuonTriggerGt20->GetBinContent(MuonTriggerGt20->FindBin(eta2, eta1));
}

float MuonAnalyzer::GetDoubleMuonTriggerSFError(pat::Muon& mu1, pat::Muon& mu2) {
    if(!isMuonTriggerFile) return 1.;
    float eta1=fabs(mu1.eta());
    float eta2=fabs(mu2.eta());
    // Muon POG enumeration is inverted 1 <-> 2
    if(mu2.pt()<20.) return MuonTriggerLt20->GetBinError(MuonTriggerLt20->FindBin(eta2, eta1));
    return MuonTriggerGt20->GetBinError(MuonTriggerGt20->FindBin(eta2, eta1));
}

float MuonAnalyzer::GetMuonIdSF(pat::Muon& mu) {
    if(!isMuonIdFile) return 1.;
    float pt=mu.pt(), eta=fabs(mu.eta());
    if(pt>=MuonPtMax) pt=MuonPtMax;
    if(pt<=MuonPtMin) pt=MuonPtMin;
    if(eta<=0.9) return MuonIdLt09->GetBinContent(MuonIdLt09->FindBin(pt));
    if(eta>0.9 && eta<=1.2) return MuonId0912->GetBinContent(MuonId0912->FindBin(pt));
    if(eta>1.2 && eta<=2.1) return MuonId1221->GetBinContent(MuonId1221->FindBin(pt));
    if(eta>2.1 && eta<=2.4) return MuonId2124->GetBinContent(MuonId2124->FindBin(pt));
    return 0.;
}

float MuonAnalyzer::GetMuonIdSFError(pat::Muon& mu) {
    if(!isMuonIdFile) return 1.;
    float pt=mu.pt(), eta=fabs(mu.eta());
    if(pt>=MuonPtMax) pt=MuonPtMax;
    if(pt<=MuonPtMin) pt=MuonPtMin;
    if(eta<=0.9) return MuonIdLt09->GetBinError(MuonIdLt09->FindBin(pt));
    if(eta>0.9 && eta<=1.2) return MuonId0912->GetBinError(MuonId0912->FindBin(pt));
    if(eta>1.2 && eta<=2.1) return MuonId1221->GetBinError(MuonId1221->FindBin(pt));
    if(eta>2.1 && eta<=2.4) return MuonId2124->GetBinError(MuonId2124->FindBin(pt));
    return 0.;
}

float MuonAnalyzer::GetMuonIsoSF(pat::Muon& mu) {
    if(!isMuonIsoFile) return 1.;
    float pt=mu.pt(), eta=fabs(mu.eta());
    if(pt>=MuonPtMax) pt=MuonPtMax;
    if(pt<=MuonPtMin) pt=MuonPtMin;
    if(eta<=0.9) return MuonIsoLt09->GetBinContent(MuonIsoLt09->FindBin(pt));
    if(eta>0.9 && eta<=1.2) return MuonIso0912->GetBinContent(MuonIso0912->FindBin(pt));
    if(eta>1.2 && eta<=2.1) return MuonIso1221->GetBinContent(MuonIso1221->FindBin(pt));
    if(eta>2.1 && eta<=2.4) return MuonIso2124->GetBinContent(MuonIso2124->FindBin(pt));
    return 0.;
}

float MuonAnalyzer::GetMuonIsoSFError(pat::Muon& mu) {
    if(!isMuonIsoFile) return 1.;
    float pt=mu.pt(), eta=fabs(mu.eta());
    if(pt>=MuonPtMax) pt=MuonPtMax;
    if(pt<=MuonPtMin) pt=MuonPtMin;
    eta=fabs(eta);
    if(eta<=0.9) return MuonIsoLt09->GetBinError(MuonIsoLt09->FindBin(pt));
    if(eta>0.9 && eta<=1.2) return MuonIso0912->GetBinError(MuonIso0912->FindBin(pt));
    if(eta>1.2 && eta<=2.1) return MuonIso1221->GetBinError(MuonIso1221->FindBin(pt));
    if(eta>2.1 && eta<=2.4) return MuonIso2124->GetBinError(MuonIso2124->FindBin(pt));
    return 0.;
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

