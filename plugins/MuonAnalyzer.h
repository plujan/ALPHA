#ifndef MUONANALYZER_H
#define MUONANALYZER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
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
        virtual float GetMuonTriggerSFIsoMu20(pat::Muon&);
        virtual float GetMuonTriggerSFErrorIsoMu20(pat::Muon&);
        virtual float GetMuonIdSFLoose(pat::Muon&);
        virtual float GetMuonIdSFLooseError(pat::Muon&);
        virtual float GetMuonIdSFTight(pat::Muon&);
        virtual float GetMuonIdSFTightError(pat::Muon&);
        virtual float GetMuonIdSFHighpt(pat::Muon&);
        virtual float GetMuonIdSFHighptError(pat::Muon&);
        virtual float GetMuonIsoSFLoose(pat::Muon&);
        virtual float GetMuonIsoSFLooseError(pat::Muon&);
        virtual float GetMuonIsoSFTight(pat::Muon&);
        virtual float GetMuonIsoSFTightError(pat::Muon&);
        virtual float GetMuonIsoSFHighpt(pat::Muon&);
        virtual float GetMuonIsoSFHighptError(pat::Muon&);
        virtual TH1F* ConvertTGraph(TGraphAsymmErrors*);
      
    private:
      
        edm::EDGetTokenT<std::vector<pat::Muon> > MuonToken;
        edm::EDGetTokenT<reco::VertexCollection> VertexToken;
	std::string MuonIdFileName;
	std::string MuonIsoFileName;
	std::string MuonHighptFileName;
	std::string MuonTriggerFileName;
	std::string DoubleMuonTriggerFileName;
        int Muon1Id, Muon2Id, Muon1Iso, Muon2Iso;
        float Muon1Pt, Muon2Pt;
        float MuonPtMax, MuonPtMin;
        
        bool isMuonTriggerFile, isDoubleMuonTriggerFile, isMuonIdFile, isMuonIsoFile, isMuonHighptFile;
        
        TFile* MuonTriggerFile;
        TFile* DoubleMuonTriggerFile;
        TFile* MuonIdFile;
        TFile* MuonHighptFile;
        TFile* MuonIsoFile;
        
        TH2F* MuonTriggerLt20;
        TH2F* MuonTriggerGt20;
	TH2F* MuonTriggerIsoMu20;
        
        TH2F* MuonIdLoose;
        TH2F* MuonIdTight;
        TH2F* MuonIdHighpt;
        TH2F* MuonIsoLoose;
        TH2F* MuonIsoTight;
        TH2F* MuonIsoHighpt;
};


#endif
