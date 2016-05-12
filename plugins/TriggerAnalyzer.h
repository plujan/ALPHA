#ifndef TRIGGERANALYZER_H
#define TRIGGERANALYZER_H

#include <iostream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
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
        TriggerAnalyzer();
        ~TriggerAnalyzer();
        
        virtual int FillTriggerBitmap(const edm::Event&, std::vector<std::string>);

      
    private:
    
};


#endif
