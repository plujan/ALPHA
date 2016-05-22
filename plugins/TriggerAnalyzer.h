#ifndef TRIGGERANALYZER_H
#define TRIGGERANALYZER_H

#include <iostream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "TFile.h"
#include "TH3.h"
#include "TKey.h"

class TriggerAnalyzer {
    public:
        TriggerAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~TriggerAnalyzer();
        
        virtual int FillTriggerBitmap(const edm::Event&, std::vector<std::string>);
      
    private:
    
        edm::EDGetTokenT<edm::TriggerResults> TriggerToken;
        std::vector<std::string> TriggerList;
};


#endif
