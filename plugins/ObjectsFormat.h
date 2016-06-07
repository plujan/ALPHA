#ifndef OBJECTSFORMAT_H
#define OBJECTSFORMAT_H

#include <string>
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
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
        static void FillPhotonType(PhotonType&, const pat::Photon*, bool);
        static void FillTauType(TauType&, const pat::Tau*, bool);
        static void FillJetType(JetType&, const pat::Jet*, bool);
        static void FillFatJetType(FatJetType&, const pat::Jet*, bool);
        static void FillMEtType(MEtType&, const pat::MET*, bool);
//        static void FillMEtCorType(MEtCorType&, const pat::MET*, bool);
        static void FillMEtFullType(MEtFullType&, const pat::MET*, bool);
        static void FillCandidateType(CandidateType&, pat::CompositeCandidate*, bool);
        static void FillLorentzType(LorentzType&, const reco::Candidate::LorentzVector*);
        static void ResetLeptonType(LeptonType&);
        static void ResetPhotonType(PhotonType&);
        static void ResetTauType(TauType&);
        static void ResetJetType(JetType&);
        static void ResetFatJetType(FatJetType&);
        static void ResetMEtType(MEtType&);
//        static void ResetMEtCorType(MEtCorType&);
        static void ResetMEtFullType(MEtFullType&);
        static void ResetCandidateType(CandidateType&);
        static void ResetLorentzType(LorentzType&);
        static std::string ListLeptonType();
        static std::string ListPhotonType();
        static std::string ListTauType();
        static std::string ListJetType();
        static std::string ListFatJetType();
        static std::string ListMEtType();
//        static std::string ListMEtCorType();
        static std::string ListMEtFullType();
        static std::string ListCandidateType();
        static std::string ListLorentzType();
        
    private:
    
};


#endif

