#ifndef JETANALYZER_H
#define JETANALYZER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "RecoParticleFlow/PFProducer/interface/Utils.h"
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

#include "RecoilCorrector.h" // From: https://github.com/cms-met/MetTools/tree/master/RecoilCorrections

#include "BTagCalibrationStandalone.h"

#include "TFile.h"
#include "TH2.h"
#include "TLorentzVector.h"

class JetAnalyzer {
    public:
        JetAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~JetAnalyzer();
        virtual std::vector<pat::Jet> FillJetVector(const edm::Event&);
        virtual void CleanJetsFromMuons(std::vector<pat::Jet>&, std::vector<pat::Muon>&, float);
        virtual void CleanJetsFromElectrons(std::vector<pat::Jet>&, std::vector<pat::Electron>&, float);
        virtual int GetNBJets(std::vector<pat::Jet>&);
        virtual pat::MET FillMetVector(const edm::Event&);
        virtual void ApplyRecoilCorrections(pat::MET&, const reco::Candidate::LorentzVector*, const reco::Candidate::LorentzVector*, int);
        virtual float GetScaleUncertainty(pat::Jet&);
        virtual float GetResolutionRatio(float);
        virtual float GetResolutionErrorUp(float);
        virtual float GetResolutionErrorDown(float);
        virtual bool isLooseJet(pat::Jet&);
        //virtual bool isMediumJet(pat::Jet&);
        virtual bool isTightJet(pat::Jet&);
        virtual bool isTightLepVetoJet(pat::Jet&);
        virtual std::vector<float> reshapeBtagDiscriminator(pat::Jet&);
      
    private:
    
        edm::EDGetTokenT<std::vector<pat::Jet> > JetToken;
        edm::EDGetTokenT<std::vector<pat::MET> > MetToken;
        edm::EDGetTokenT<reco::JetCorrector> CorToken;
        edm::EDGetTokenT<edm::ValueMap<float>> QGToken;
        int JetId;
        float Jet1Pt, Jet2Pt, JetEta;
        bool AddQG, RecalibrateJets, RecalibrateMass;
        std::string BTag;
        int Jet1BTag, Jet2BTag;
        std::string BTagDB;
        bool UseRecoil;
        std::string RecoilMCFile;
        std::string RecoilDataFile;
        
        bool isJESFile;
        
        TFile* JESFile;
        TH2F* hist;
        
        RecoilCorrector* recoilCorr;
};

#endif
