#ifndef MUONANALYZER_H
#define MUONANALYZER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TFile.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"

class MuonAnalyzer {
    public:
        MuonAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~MuonAnalyzer();
        virtual std::vector<pat::Muon> FillMuonVector(const edm::Event&);
        virtual bool IsCustomTracker(pat::Muon&, const reco::Vertex*);
        virtual float GetDoubleMuonTriggerSF(pat::Muon&, pat::Muon&);
        virtual float GetDoubleMuonTriggerSFError(pat::Muon&, pat::Muon&);
        virtual float GetMuonIdSF(pat::Muon&);
        virtual float GetMuonIdSFError(pat::Muon&);
        virtual float GetMuonIsoSF(pat::Muon&);
        virtual float GetMuonIsoSFError(pat::Muon&);
        virtual TH1F* ConvertTGraph(TGraphAsymmErrors*);
      
    private:
      
        edm::EDGetTokenT<std::vector<pat::Muon> > MuonToken;
        edm::EDGetTokenT<reco::VertexCollection> VertexToken;
        int Muon1Id, Muon2Id, Muon1Iso, Muon2Iso;
        float Muon1Pt, Muon2Pt;
        float MuonPtMax, MuonPtMin;
        
        bool isMuonTriggerFile, isMuonIdFile, isMuonIsoFile;
        
        TFile* MuonTriggerFile;
        TFile* MuonIdFile;
        TFile* MuonIsoFile;
        
        TH2F* MuonTriggerLt20;
        TH2F* MuonTriggerGt20;
        
        TH1F* MuonIdLt09;
        TH1F* MuonId0912;
        TH1F* MuonId1221;
        TH1F* MuonId2124;
        TH1F* MuonIsoLt09;
        TH1F* MuonIso0912;
        TH1F* MuonIso1221;
        TH1F* MuonIso2124;
};


#endif
