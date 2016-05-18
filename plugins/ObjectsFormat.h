#ifndef OBJECTSFORMAT_H
#define OBJECTSFORMAT_H

#include <string>
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "Objects.h"

class ObjectsFormat {
    
    public:
        ObjectsFormat() {};
        ~ObjectsFormat() {};
        
        static void FillElectronType(LeptonType&, const pat::Electron*, bool);
        static void FillMuonType(LeptonType&, const pat::Muon*, bool);
        static void FillJetType(JetType&, const pat::Jet*, bool);
        static void FillMEtType(MEtType&, const pat::MET*, bool);
        static void FillMEtCorType(MEtCorType&, const pat::MET*, bool);
        static void FillMEtFullType(MEtFullType&, const pat::MET*, bool);
        static void ResetLeptonType(LeptonType&);
        static void ResetJetType(JetType&);
        static void resetMEtType(MEtType&);
        static std::string ListLeptonType();
        static std::string ListJetType();
        
    private:
    
};


#endif

