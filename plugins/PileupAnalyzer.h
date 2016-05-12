#ifndef PILEUPANALYZER_H
#define PILEUPANALYZER_H

#include <iostream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TFile.h"
#include "TH1.h"

class PileupAnalyzer {
    public:
        PileupAnalyzer();
        ~PileupAnalyzer();
        virtual float GetPUWeight(const edm::Event&);

      
    private:

        edm::LumiReWeighting* LumiWeights;
};

#endif
