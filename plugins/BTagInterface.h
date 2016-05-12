#ifndef BTAGINTERFACE_H
#define BTAGINTERFACE_H

#include <iostream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "TFile.h"
#include "TH2.h"
#include "TGraph.h"

class BTagInterface {
  public:
    BTagInterface(std::string);
    ~BTagInterface();
    virtual void FillBTagVector(const edm::Event&, std::vector<pat::Jet>&);
    virtual float GetReshapedDiscr(pat::Jet&, int=0);
    virtual float GetNewWorkingPoint(int, int, float, float, int);
    virtual float GetWorkingPoint(int);
    virtual float GetScaleFactor(int, float);
    virtual float GetScaleFactor(int, float, float, int=0);
    virtual float GetScaleFactorMistag(int, float, float, int=0);
    
  private:
    std::string Tagger;
    float WorkingPoint[5];
    float Integral[3][5];
    bool isBDFile, isShapesFile;
    
    TFile* BDFile;
    TFile* ShapesFile;
    TH2F* SFb[3];
    TH2F* SFl[3];
    TH1F* Shape[3]; // 0: B, 1: C, 2: L
};

#endif
