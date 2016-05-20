#ifndef ELECTRONANALYZER_H
#define ELECTRONANALYZER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "TFile.h"
#include "TH2.h"

class ElectronAnalyzer {
    public:
        ElectronAnalyzer(const edm::ParameterSet&, edm::ConsumesCollector&&);
        ~ElectronAnalyzer();
        virtual std::vector<pat::Electron> FillElectronVector(const edm::Event&);
        //virtual bool IsElectronCut(pat::Electron&, const reco::Vertex*, int);
        //virtual bool IsLooseElectron(pat::Electron&, const reco::Vertex*);
        //virtual bool IsLooseMVAElectron(pat::Electron&);
        //virtual bool IsTightMVAElectron(pat::Electron&);
        virtual float GetDoubleElectronTriggerSF(pat::Electron&, pat::Electron&);
        virtual float GetElectronIdSF(pat::Electron&);
        virtual float GetElectronIdSFError(pat::Electron&);
        virtual float GetElectronIsoSF(pat::Electron&);
        virtual float GetElectronIsoSFError(pat::Electron&);
        virtual float GetLooseElectronSF(pat::Electron&);
      
    private:
      
        edm::EDGetTokenT<std::vector<pat::Electron> > ElectronToken;
        edm::EDGetTokenT<reco::VertexCollection> VertexToken;
	edm::EDGetTokenT<edm::ValueMap<bool>> EleVetoIdMapToken;
	edm::EDGetTokenT<edm::ValueMap<bool>> EleLooseIdMapToken;
	edm::EDGetTokenT<edm::ValueMap<bool>> EleMediumIdMapToken;
	edm::EDGetTokenT<edm::ValueMap<bool>> EleTightIdMapToken;
	edm::EDGetTokenT<edm::ValueMap<bool>> EleHEEPIdMapToken;
	edm::EDGetTokenT<edm::ValueMap<bool>> EleMVANonTrigMediumIdMapToken;
	edm::EDGetTokenT<edm::ValueMap<bool>> EleMVANonTrigTightIdMapToken;
	edm::EDGetTokenT<edm::ValueMap<bool>> EleMVATrigMediumIdMapToken;
	edm::EDGetTokenT<edm::ValueMap<bool>> EleMVATrigTightIdMapToken;
        int Electron1Id, Electron2Id, Electron1Iso, Electron2Iso;
        float Electron1Pt, Electron2Pt;
        float EleTriggerPtMax, EleEfficiencyPtMax;
        
        bool isEleTriggerFile, isEleEfficiencyFile;
        
        TFile* EleTriggerFile;
        TFile* EleEfficiencyFile;
        
        TH2F* EleTriggerDATAHighLeg;
        TH2F* EleTriggerDATALowLeg;
        TH2F* EleTriggerMCHighLeg;
        TH2F* EleTriggerMCLowLeg;
        TH2F* ElectronId;
        TH2F* ElectronIso;
};


#endif
