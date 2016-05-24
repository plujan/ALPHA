#ifndef PHOTONANALYZER_H
#define PHOTONANALYZER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TFile.h"
#include "TH2.h"

class PhotonAnalyzer {
    public:
        PhotonAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~PhotonAnalyzer();
        virtual std::vector<pat::Photon> FillPhotonVector(const edm::Event&);
        //virtual bool isLoosePhoton(pat::Photon&, const reco::Vertex*);
      
    private:
      
        edm::EDGetTokenT<std::vector<pat::Photon> > PhotonToken;
        edm::EDGetTokenT<reco::VertexCollection> VertexToken;
	//edm::EDGetTokenT<edm::ValueMap<bool>> PhoLooseIdMapToken;
	//edm::EDGetTokenT<edm::ValueMap<bool>> PhoMediumIdMapToken;
	//edm::EDGetTokenT<edm::ValueMap<bool>> PhoTightIdMapToken;
	//edm::EDGetTokenT<edm::ValueMap<bool>> PhoMVANonTrigMediumIdMapToken;
        int Photon1Id, Photon2Id;//, Photon1Iso, Photon2Iso;
        float Photon1Pt, Photon2Pt;
    
};


#endif
