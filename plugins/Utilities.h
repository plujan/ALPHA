#ifndef UTILITIES_H
#define UTILITIES_H


// dR and dPhi
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefToBase.h"

// Jets and PF
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "RecoParticleFlow/PFProducer/interface/Utils.h"
// Transient Track and IP
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
// Gen Info
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
// Association
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

// Pat
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
// Trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
// KinFitter
#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TFitParticlePtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEScaledMomDev.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"

#include "Numbers.h"

#include "TLorentzVector.h"
#include "TVector3.h"


//bool SortByPt(pat::Jet j1, pat::Jet j2) {return(j1.pt() > j2.pt());}
//bool SortByCSV(pat::Jet j1, pat::Jet j2) {return(j1.bDiscriminator("combinedSecondaryVertexBJetTags") > j2.bDiscriminator("combinedSecondaryVertexBJetTags"));}
//bool SortByCSVV1(pat::Jet j1, pat::Jet j2) {return(j1.bDiscriminator("combinedSecondaryVertexV1BJetTags") > j2.bDiscriminator("combinedSecondaryVertexV1BJetTags"));}
//bool SortByCSVSL(pat::Jet j1, pat::Jet j2) {return(j1.bDiscriminator("combinedSecondaryVertexSoftPFLeptonV1BJetTags") > j2.bDiscriminator("combinedSecondaryVertexSoftPFLeptonV1BJetTags"));}
//bool SortByCSVV1Reshaped(pat::Jet j1, pat::Jet j2) {return(j1.bDiscriminator("combinedSecondaryVertexV1BJetTagsReshaped") > j2.bDiscriminator("combinedSecondaryVertexV1BJetTagsReshaped"));}

class Utilities {
  public:
    Utilities() {};
    ~Utilities() {};
    
//    static bool isElectronCut(const pat::Electron*, const reco::Vertex*, float, int);
//    static bool isLooseElectron(const pat::Electron*, const reco::Vertex*, float);
//    static bool isHZZElectron(const pat::Electron*);
//    static bool isHWWElectron(const pat::Electron*);
//    
//    static bool isLooseMuon(const pat::Muon*);
//    static bool isSoftMuon (const pat::Muon*, const reco::Vertex*);
//    static bool isTightMuon(const pat::Muon*, const reco::Vertex*);
//    
//    static bool isLooseJet (const pat::Jet*);
//    static bool isMediumJet(const pat::Jet*);
//    static bool isTightJet (const pat::Jet*);
//    static void addLHEInfoToJet(pat::Jet&, std::vector<reco::Candidate::LorentzVector>&);
//    static int GetJetFlavour(const pat::Jet*);
//    
//    static int  FindMotherId(const reco::Candidate*);
//    // Event Shape
//    static void EventShape(std::vector<pat::Jet>*, float&, float&);
//    static std::vector<float> ReturnEventShape(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, bool);
//    static float ReturnCentrality(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
//    // Radiation
//    static void ReturnThetaPull(const pat::Jet*, const pat::Jet*, float&, float&, float&, float&);
//    static float ReturnGirth(const pat::Jet*);
    // Angular
    static float ReturnCosThetaStar(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
    static float ReturnCosTheta1(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
    static float ReturnCosTheta2(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
    static float ReturnPhi(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
    static float ReturnPhi1(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&);
    // Kinematics and reconstruction
    static float PerformKinematicFit(pat::Jet*, pat::Jet*, reco::Candidate::LorentzVector*, reco::Candidate::LorentzVector*, float);
    static float RecoverNeutrinoPz(const reco::Particle::LorentzVector*, const reco::Particle::LorentzVector*);
    // Miscellanea
    static pat::CompositeCandidate RecoilMassFormula(pat::CompositeCandidate&, pat::MET&);
    static int FindMotherId(const reco::GenParticle*);
    
    
  private:
    
};

/*
// --------------------------------------------------
// -------------------- ELECTRONS -------------------
// --------------------------------------------------

// Electron Cut-based Quality ID: see https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification#Electron_ID_Working_Points
bool Utilities::isElectronCut(const pat::Electron* el, const reco::Vertex* vertex, float iso, int q) {
  // Electron Quality: 0 - Veto, 1 - Loose, 2 - Medium, 3 - Tight
  //                                VBTF90     VBTF 80    VBTF 70
  if(q<0 || q>3) return false;
  //     ECAL BARREL                                   ECAL ENDCAP
  float sieieEB[4]    ={0.01,  0.01,  0.01,  0.01 },  sieieEE[4]    ={0.03,  0.03,  0.03,  0.03 };
  float dEtaInEB[4]   ={0.007, 0.007, 0.004, 0.004},  dEtaInEE[4]   ={0.010, 0.009, 0.007, 0.005};
  float dPhiInEB[4]   ={0.8,   0.15,  0.06,  0.03 },  dPhiInEE[4]   ={0.7,   0.10,  0.03,  0.02 };
  float hOeEB[4]      ={0.15,  0.12,  0.12,  0.12 },  hOeEE[4]      ={1.e10, 0.10,  0.10,  0.10 };
  float d0vtxEB[4]    ={0.04,  0.02,  0.02,  0.02 },  d0vtxEE[4]    ={0.04,  0.02,  0.02,  0.02 };
  float dZvtxEB[4]    ={0.2,   0.2,   0.1,   0.1  },  dZvtxEE[4]    ={0.2,   0.2,   0.1,   0.1  };
  float abs1oE1opEB[4]={0.,    0.05,  0.05,  0.05 },  abs1oE1opEE[4]={0.,    0.05,  0.05,  0.05 };
//float pfIsoOPtEB[4] ={0.15,  0.15,  0.15,  0.10 },  pfIsoOPtEE[4] ={0.15,  0.15,  0.15,  0.10 };
//float vtxFitEB[4]   ={1.e-6, 1.e-6, 1.e-6, 1.e-6},  vtxFitEE[4]   ={1.e-6, 1.e-6, 1.e-6, 1.e-6};
  float mHitsEB[4]    ={1.,    1.,    1.,    0.   },  mHitsEE[4]    ={1.,    1.,    1.,    0.   };
  
  bool isEB=el->isEB();
  bool isEE=el->isEE();
  
  float sieie=el->scSigmaIEtaIEta();
  float dEtaIn=fabs(el->deltaEtaSuperClusterTrackAtVtx());
  float dPhiIn=fabs(el->deltaPhiSuperClusterTrackAtVtx());
  float hOe=el->hadronicOverEm();
  float d0vtx=fabs( el->gsfTrack()->d0() - vertex->x()*sin(el->gsfTrack()->phi()) + vertex->y()*cos(el->gsfTrack()->phi()) );
  float dZvtx=fabs( (el->vz()-vertex->z()) -( (el->vx()-vertex->x())*el->px()+(el->vy()-vertex->y())*el->py() )/el->pt()*el->pz()/el->pt() );
  float E = el->superCluster()->energy();
  float p = 1./(el->eSuperClusterOverP()/E);
  float abs1oE1op=fabs(1./E - 1./p);
  float mHits=el->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
  
  if(isEB) if(sieie>sieieEB[q] || dEtaIn>dEtaInEB[q] || dPhiIn>dPhiInEB[q] || hOe>hOeEB[q] || d0vtx>d0vtxEB[q] || dZvtx>dZvtxEB[q] || abs1oE1op>abs1oE1opEB[q] || mHits>mHitsEB[q]) return false;
  if(isEE) if(sieie>sieieEE[q] || dEtaIn>dEtaInEE[q] || dPhiIn>dPhiInEE[q] || hOe>hOeEE[q] || d0vtx>d0vtxEE[q] || dZvtx>dZvtxEE[q] || abs1oE1op>abs1oE1opEE[q] || mHits>mHitsEE[q]) return false;
  return true;
}

// https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification#Glossary_of_Variables
bool Utilities::isLooseElectron(const pat::Electron* el, const reco::Vertex* vertex, float iso) {
  if(el->isEB()) {
    if(iso<0.15
      && fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.007
      && fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.15
      && el->scSigmaIEtaIEta()<0.01
      && el->hadronicOverEm()<0.12
      && el->gsfTrack()->dxy(vertex->position())<0.02
      && el->gsfTrack()->dz(vertex->position())<0.2
      && fabs(1./el->ecalEnergy() - 1./el->trackMomentumAtVtx().R())<0.05
//      && el->reco::GsfElectron::fbrem()<10.e-6
      && el->gsfTrack()->trackerExpectedHitsInner().numberOfHits()<=1
    ) return true;
  }
  else if(el->isEE()) {
    if(iso<0.15
      && fabs(el->deltaEtaSuperClusterTrackAtVtx())<0.009
      && fabs(el->deltaPhiSuperClusterTrackAtVtx())<0.10
      && el->scSigmaIEtaIEta()<0.03
      && el->hadronicOverEm()<0.10
      && el->gsfTrack()->dxy(vertex->position())<0.02
      && el->gsfTrack()->dz(vertex->position())<0.2
      && fabs(1./el->ecalEnergy() - 1./el->trackMomentumAtVtx().R())<0.05
//      && el->reco::GsfElectron::fbrem()<10.e-6
      && el->gsfTrack()->trackerExpectedHitsInner().numberOfHits()<=1
    ) return true;
  }
  return false;
}


bool Utilities::isHZZElectron(const pat::Electron* el) {
  //if(el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > 0) return false;
  float mva=el->electronID("mvaTrigV0");
  float eta=fabs(el->eta());
  // AN2013_108_v8 page 31
  if(el->pt()>=10.) {
    if(eta<0.8)               {if(mva<-0.34) return false;}
    if(eta>=0.8 && eta<1.479) {if(mva<-0.65) return false;}
    if(eta>=1.479)            {if(mva<+0.60) return false;}
  }
  else {
    if(eta<0.8)               {if(mva<+0.47) return false;}
    if(eta>=0.8 && eta<1.479) {if(mva<0.004) return false;}
    if(eta>=1.479)            {if(mva<0.295) return false;}
  }
  // https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification#Recommended_Working_Points_With
  //if(eta<0.8)               {if(mva<0.94) return false;}
  //if(eta>=0.8 && eta<1.479) {if(mva<0.85) return false;}
  //if(eta>=1.479)            {if(mva<0.92) return false;}
  //  Old MVA cuts https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMultivariateIsoElectrons
  //if(eta<0.8)               {if(mva<0.735) return false;}
  //if(eta>=0.8 && eta<1.479) {if(mva<0.467) return false;}
  //if(eta>=1.479)            {if(mva<0.795) return false;}
  // VH
  //return ((eta<0.8 && mva>0.913 && iso<0.105) || (eta>=0.8 && eta<1.479 && mva>0.964 && iso<0.178) || (eta>=1.479 && eta<2.5 && mva>0.899 && iso<0.150)); // WP80
  return true;
}

// https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification#Recommended_Working_Points_With
bool Utilities::isHWWElectron(const pat::Electron* el) {
  // Conversion Veto
  if(!el->passConversionVeto()) return false;
  // Missing hits veto
  if(el->gsfTrack()->trackerExpectedHitsInner().numberOfHits()>0) return false;
  // Triggering MVA Cut
  float mva=el->electronID("mvaTrigV0");
  float eta=fabs(el->eta());
  if(el->pt()<20.) {
    if(eta<0.8)               {if(mva<+0.00) return false;}
    if(eta>=0.8 && eta<1.479) {if(mva<+0.10) return false;}
    if(eta>=1.479)            {if(mva<+0.62) return false;}
  }
  else {
    if(eta<0.8)               {if(mva<+0.94) return false;}
    if(eta>=0.8 && eta<1.479) {if(mva<+0.85) return false;}
    if(eta>=1.479)            {if(mva<+0.92) return false;}
  }
  return true;
}
// --------------------------------------------------
// --------------------   MUONS   -------------------
// --------------------------------------------------

// Muon Quality ID: see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
bool Utilities::isLooseMuon(const pat::Muon* muon) {
  return muon->isPFMuon() && (muon->isGlobalMuon() || muon->isTrackerMuon());
}
bool Utilities::isSoftMuon(const pat::Muon* muon, const reco::Vertex* vertex) {
  return  muon::isGoodMuon(*muon, muon::TMOneStationTight)
    && muon->track()->hitPattern().trackerLayersWithMeasurement()       > 5
    && muon->innerTrack()->hitPattern().pixelLayersWithMeasurement()    > 1
    && muon->muonBestTrack()->normalizedChi2()                          < 1.8
    && muon->innerTrack()->dxy(vertex->position())                      < 3.
    && muon->innerTrack()->dz(vertex->position())                       < 30.
  ;
}
bool Utilities::isTightMuon(const pat::Muon* muon, const reco::Vertex* vertex) {
  return  muon->isGlobalMuon()
    && muon->isPFMuon()
    && muon->muonBestTrack()->normalizedChi2()                          < 10.
    && (muon->globalTrack().isNonnull() ? muon->globalTrack()->hitPattern().numberOfValidMuonHits() : -1)        > 0
    && muon->numberOfMatchedStations()                                  > 1
    && fabs(muon->muonBestTrack()->dxy(vertex->position()))             < 0.2 
    && fabs(muon->muonBestTrack()->dz(vertex->position()))              < 0.5
    && muon->innerTrack()->hitPattern().numberOfValidPixelHits()        > 0
    && muon->track()->hitPattern().trackerLayersWithMeasurement()       > 5
  ;
}

// --------------------------------------------------
// --------------------   JETS    -------------------
// --------------------------------------------------

// PFJet Quality ID: see https://twiki.cern.ch/twiki/bin/view/CMS/JetID
bool Utilities::isLooseJet(const pat::Jet* jet) {
  if(jet->neutralHadronEnergyFraction()>0.99) return false;
  if(jet->neutralEmEnergyFraction()>0.99) return false;
  if(jet->numberOfDaughters()<=1) return false;
  if(fabs(jet->eta())<2.4) {
    if(jet->chargedHadronEnergyFraction()<=0.) return false;
    if(jet->chargedEmEnergyFraction()>0.99) return false;
    if(jet->chargedMultiplicity()<=0) return false;
  }
  return true;
}
bool Utilities::isMediumJet(const pat::Jet* jet) {
  if(jet->neutralHadronEnergyFraction()>0.95) return false;
  if(jet->neutralEmEnergyFraction()>0.95) return false;
  if(jet->numberOfDaughters()<=1) return false;
  if(fabs(jet->eta())<2.4) {
    if(jet->chargedHadronEnergyFraction()<=0.) return false;
    if(jet->chargedEmEnergyFraction()>0.99) return false;
    if(jet->chargedMultiplicity()<=0) return false;
  }
  return true;
}
bool Utilities::isTightJet(const pat::Jet* jet) {
  if(jet->neutralHadronEnergyFraction()>0.90) return false;
  if(jet->neutralEmEnergyFraction()>0.90) return false;
  if(jet->numberOfDaughters()<=1) return false;
  if(fabs(jet->eta())<2.4) {
    if(jet->chargedHadronEnergyFraction()<=0.) return false;
    if(jet->chargedEmEnergyFraction()>0.99) return false;
    if(jet->chargedMultiplicity()<=0) return false;
  }
  return true;
}

void Utilities::addLHEInfoToJet(pat::Jet& jet, std::vector<reco::Candidate::LorentzVector>& lhe) {
  if(!jet.genParton()) return;
  float dRTmp(99.), dRMin(99.);
  for(unsigned int i=0; i<lhe.size(); i++) {
    dRTmp=deltaR(lhe.at(i), *jet.genParton());
    if(dRTmp<dRMin) {
      dRMin=dRTmp;
      jet.addUserFloat("ptLhe", (float)lhe.at(i).pt());
      jet.addUserFloat("etaLhe", (float)lhe.at(i).eta());
      jet.addUserFloat("phiLhe", (float)lhe.at(i).phi());
    }
  }
}

int Utilities::GetJetFlavour(const pat::Jet* jet) {
  if(abs(jet->partonFlavour())==5) return 0;
  if(abs(jet->partonFlavour())==4) return 1;
  if(abs(jet->partonFlavour())<4 || abs(jet->partonFlavour())==21) return 2;
  return -1;
}



// --------------------------------------------------
// --------------------   OTHERS  -------------------
// --------------------------------------------------

int Utilities::FindMotherId(const reco::Candidate* p) {
  if(p->numberOfMothers()<1) return 0;
  int dId, mId;
  do {
    dId=abs(p->pdgId());
    try {
      p=p->mother();
      mId=abs(p->pdgId());
    }
    catch(...) {
      return -1;
    }
  } while(p->numberOfMothers()>0 && dId==mId);
  return abs(p->pdgId());
}


// Event Shape Variables, following http://home.fnal.gov/~mrenna/lutp0613man2/node233.html
void Utilities::EventShape(std::vector<pat::Jet>* Jets, float& sphericity, float& thrust) {
  TLorentzVector Ptot;
  for(std::vector<pat::Jet>::const_iterator jet=Jets->begin(); jet!=Jets->end(); ++jet) {
    Ptot(0)+=jet->px();
    Ptot(1)+=jet->py();
    Ptot(2)+=jet->pz();
    Ptot(3)+=jet->energy();
  }

  TVector3 beta=Ptot.BoostVector();

  TMatrixD PTensor(3,3);
  double p2=0.0;
  for(std::vector<pat::Jet>::const_iterator jet=Jets->begin(); jet!=Jets->end(); ++jet) {
    TLorentzVector pj(jet->px(),jet->py(),jet->pz(),jet->p());
    pj.Boost(-beta);
    p2+=pj.P()*pj.P();
    PTensor(0,0)+=pj.Px()*pj.Px();
    PTensor(0,1)+=pj.Px()*pj.Py();
    PTensor(0,2)+=pj.Px()*pj.Pz();
    PTensor(1,0)+=pj.Py()*pj.Px();
    PTensor(1,1)+=pj.Py()*pj.Py();
    PTensor(1,2)+=pj.Py()*pj.Pz();
    PTensor(2,0)+=pj.Pz()*pj.Px();
    PTensor(2,1)+=pj.Pz()*pj.Py();
    PTensor(2,2)+=pj.Pz()*pj.Pz();
  }
  PTensor(0,0)/=p2;
  PTensor(0,1)/=p2;
  PTensor(0,2)/=p2;
  PTensor(1,0)/=p2;
  PTensor(1,1)/=p2;
  PTensor(1,2)/=p2;
  PTensor(2,0)/=p2;
  PTensor(2,1)/=p2;
  PTensor(2,2)/=p2;

  TVectorD EigenVal(3);
  TMatrixD EigenVec(3,3);
  EigenVec=PTensor.EigenVectors(EigenVal);

  sphericity=3./2.*(EigenVal[1]+EigenVal[2]);
  thrust=3./2.*EigenVal[2];
}



// Event Shape Variables, following http://home.fnal.gov/~mrenna/lutp0613man2/node233.html
std::vector<float> Utilities::ReturnEventShape(const reco::Candidate::LorentzVector& L1, const reco::Candidate::LorentzVector& L2, const reco::Candidate::LorentzVector& B1, const reco::Candidate::LorentzVector& B2, bool ferox) {
  std::vector<TLorentzVector> P4;
  P4.push_back(TLorentzVector(L1.px(), L1.py(), L1.pz(), L1.energy()));
  P4.push_back(TLorentzVector(L2.px(), L2.py(), L2.pz(), L2.energy()));
  P4.push_back(TLorentzVector(B1.px(), B1.py(), B1.pz(), B1.energy()));
  P4.push_back(TLorentzVector(B2.px(), B2.py(), B2.pz(), B2.energy()));
  
  TLorentzVector P4tot;
  for(unsigned i=0; i<P4.size(); i++) P4tot+=P4.at(i);
  for(unsigned i=0; i<P4.size(); i++) P4tot.Boost(-P4tot.BoostVector());
  
  float den(0.);
  TMatrixF STensor(3, 3);
  for(unsigned i=0; i<P4.size(); i++) {
    if(ferox) den+=P4.at(i).M();
    else den+=P4.at(i).M2();
    STensor(0, 0)+=P4.at(i).Px()*P4.at(i).Px();
    STensor(0, 1)+=P4.at(i).Px()*P4.at(i).Py();
    STensor(0, 2)+=P4.at(i).Px()*P4.at(i).Pz();
    STensor(1, 0)+=P4.at(i).Py()*P4.at(i).Px();
    STensor(1, 1)+=P4.at(i).Py()*P4.at(i).Py();
    STensor(1, 2)+=P4.at(i).Py()*P4.at(i).Pz();
    STensor(2, 0)+=P4.at(i).Pz()*P4.at(i).Px();
    STensor(2, 1)+=P4.at(i).Pz()*P4.at(i).Py();
    STensor(2, 2)+=P4.at(i).Pz()*P4.at(i).Pz();
    if(ferox) for(int m=0; m<3; m++) for(int n=0; n<3; n++) STensor(m, n)/=P4.at(i).M();
  }
  for(int m=0; m<3; m++) for(int n=0; n<3; n++) STensor(m, n)/=den;

  TVectorF EigenVal(3);
  TMatrixF EigenVec(3,3);
  EigenVec=STensor.EigenVectors(EigenVal);
  float Norm(EigenVal[0]+EigenVal[1]+EigenVal[2]);
  std::vector<float> EV(3);
  for(int i=0; i<3; i++) EV[i]=EigenVal[i]/Norm;
  std::sort(EV.begin(), EV.end());
  return EV;// Warning, Eigenvalues ordered decreasingly istead of increasingly as in paper http://arxiv.org/abs/1010.3698
}

float Utilities::ReturnCentrality(const reco::Candidate::LorentzVector& L1, const reco::Candidate::LorentzVector& L2, const reco::Candidate::LorentzVector& B1, const reco::Candidate::LorentzVector& B2) {
  TLorentzVector tL1(L1.px(), L1.py(), L1.pz(), L1.energy());
  TLorentzVector tL2(L2.px(), L2.py(), L2.pz(), L2.energy());
  TLorentzVector tB1(B1.px(), B1.py(), B1.pz(), B1.energy());
  TLorentzVector tB2(B2.px(), B2.py(), B2.pz(), B2.energy());
  TLorentzVector tA(tL1+tL2+tB1+tB2);
  // Boost objects to the A rest frame
  tL1.Boost( -tA.BoostVector() );
  tL2.Boost( -tA.BoostVector() );
  tB1.Boost( -tA.BoostVector() );
  tB2.Boost( -tA.BoostVector() );
  // Calc Centrality
  float value = (tL1.Pt()+tL2.Pt()+tB1.Pt()+tB2.Pt())/(tL1.Energy()+tL2.Energy()+tB1.Energy()+tB2.Energy());
  if(value!=value || isinf(value) || value<0.) return 0.;
  return value;
}

// --------------------   RADIATION  -------------------

// Pull angle http://arxiv.org/abs/1001.5027v3
void Utilities::ReturnThetaPull(const pat::Jet* Jet1, const pat::Jet* Jet2, float& s_pull1, float& s_pull2, float& b_pull1, float& b_pull2) {
  // Jet direction only with charged candidates
  reco::Candidate::LorentzVector Axis1;
  for(unsigned int i=0; i<Jet1->getPFConstituents().size(); i++) {
    if(Jet1->getPFConstituent(i)->charge()!=0) {
      Axis1+=Jet1->getPFConstituent(i)->p4();
    }
  }
  reco::Candidate::LorentzVector Axis2;
  for(unsigned int i=0; i<Jet2->getPFConstituents().size(); i++) {
    if(Jet2->getPFConstituent(i)->charge()!=0) {
      Axis2+=Jet2->getPFConstituent(i)->p4();
    }
  }
  // t Vector, lies int the y-Phi plane, it is not a 4-vector: see paper
  TVector2 t1(0, 0), t2(0, 0);
  for(unsigned int i=0; i<Jet1->getPFConstituents().size(); i++) {
    if(Jet1->getPFConstituent(i)->charge()!=0) {
      TVector2 r(Jet1->getPFConstituent(i)->rapidity()-Axis1.Rapidity(), deltaPhi(Jet1->getPFConstituent(i)->phi(), Axis1.Phi()));
      t1+=( Jet1->getPFConstituent(i)->pt() * r.Mod() / Axis1.pt() ) * r;
    }
  }
  for(unsigned int i=0; i<Jet2->getPFConstituents().size(); i++) {
    if(Jet2->getPFConstituent(i)->charge()!=0) {
      TVector2 r(Jet2->getPFConstituent(i)->rapidity()-Axis2.Rapidity(), deltaPhi(Jet2->getPFConstituent(i)->phi(), Axis2.Phi()));
      t2+=( Jet2->getPFConstituent(i)->pt() * r.Mod() / Axis2.pt() ) * r;
    }
  }
  // Signal
  // Axis12: Jet1->Jet2
  TVector2 Axis12( Axis2.Rapidity()-Axis1.Rapidity(), deltaPhi(Axis2.Phi(), Axis1.Phi()) );
  if(t1.Mod()>0.) s_pull1 = t1.DeltaPhi( Axis12 );
  if(t2.Mod()>0.) s_pull2 = t2.DeltaPhi(-1*Axis12);
  
  // Background
  // Beams
  TVector2 Beam1( Axis1.Rapidity()>Axis2.Rapidity() ? +1. : -1 , 0.);
  TVector2 Beam2( -1*Beam1 );
  if(t1.Mod()>0.) b_pull1 = t1.DeltaPhi( Beam1 );
  if(t2.Mod()>0.) b_pull2 = t2.DeltaPhi( Beam2 );
}

float Utilities::ReturnGirth(const pat::Jet* Jet) {
  float girth(0.);
  reco::Candidate::LorentzVector Axis;
  for(unsigned int i=0; i<Jet->getPFConstituents().size(); i++) {
    if(Jet->getPFConstituent(i)->charge()!=0) {
      Axis+=Jet->getPFConstituent(i)->p4();
    }
  }
  for(unsigned int i=0; i<Jet->getPFConstituents().size(); i++) {
    if(Jet->getPFConstituent(i)->charge()!=0) {
      TVector2 r(Jet->getPFConstituent(i)->rapidity()-Axis.Rapidity(), deltaPhi(Jet->getPFConstituent(i)->phi(), Axis.Phi()));
      girth+=Jet->getPFConstituent(i)->pt() * r.Mod() / Axis.pt();
    }
  }
  return girth;
}
*/

// --------------------   ANGULAR  -------------------

float Utilities::ReturnCosThetaStar(const reco::Candidate::LorentzVector& theA, const reco::Candidate::LorentzVector& theZ) {
  TLorentzVector pA(theA.px(), theA.py(), theA.pz(), theA.energy());
  TLorentzVector pZ(theZ.px(), theZ.py(), theZ.pz(), theZ.energy());
  // Boost the Z to the A rest frame
  pZ.Boost( -pA.BoostVector() );
  float value=pZ.CosTheta();
  if(value!=value || isinf(value)) return 0.;
  return value;
}

float Utilities::ReturnCosTheta1(const reco::Candidate::LorentzVector& theA, const reco::Candidate::LorentzVector& theL1, const reco::Candidate::LorentzVector& theL2, const reco::Candidate::LorentzVector& theB1, const reco::Candidate::LorentzVector& theB2) {
  TLorentzVector pA(theA.px(), theA.py(), theA.pz(), theA.energy());
  TLorentzVector pL1(theL1.px(), theL1.py(), theL1.pz(), theL1.energy());
  TLorentzVector pL2(theL2.px(), theL2.py(), theL2.pz(), theL2.energy());
  TLorentzVector pB1(theB1.px(), theB1.py(), theB1.pz(), theB1.energy());
  TLorentzVector pB2(theB2.px(), theB2.py(), theB2.pz(), theB2.energy());
  // Boost objects to the A rest frame
  pL1.Boost( -pA.BoostVector() );
  pL2.Boost( -pA.BoostVector() );
  pB1.Boost( -pA.BoostVector() );
  pB2.Boost( -pA.BoostVector() );
  // Reconstruct H in A rest frame
  TLorentzVector pHr = pB1 + pB2;
  // cos theta = H dot L1 / (|H|*|L1|)
  float value=pHr.Vect().Dot( pL1.Vect() ) / ( pHr.Vect().Mag()*pL1.Vect().Mag() );
  if(value!=value || isinf(value)) return 0.;
  return value;
}

float Utilities::ReturnCosTheta2(const reco::Candidate::LorentzVector& theA, const reco::Candidate::LorentzVector& theL1, const reco::Candidate::LorentzVector& theL2, const reco::Candidate::LorentzVector& theB1, const reco::Candidate::LorentzVector& theB2) {
  TLorentzVector pA(theA.px(), theA.py(), theA.pz(), theA.energy());
  TLorentzVector pL1(theL1.px(), theL1.py(), theL1.pz(), theL1.energy());
  TLorentzVector pL2(theL2.px(), theL2.py(), theL2.pz(), theL2.energy());
  TLorentzVector pB1(theB1.px(), theB1.py(), theB1.pz(), theB1.energy());
  TLorentzVector pB2(theB2.px(), theB2.py(), theB2.pz(), theB2.energy());
  // Boost objects to the A rest frame
  pL1.Boost( -pA.BoostVector() );
  pL2.Boost( -pA.BoostVector() );
  pB1.Boost( -pA.BoostVector() );
  pB2.Boost( -pA.BoostVector() );
  // Reconstruct Z in A rest frame
  TLorentzVector pZr = pL1 + pL2;
  // cos theta = Z dot B1 / (|Z|*|B1|)
  float value=pZr.Vect().Dot( pB1.Vect() ) / ( pZr.Vect().Mag()*pB1.Vect().Mag() );
  if(value!=value || isinf(value)) return 0.;
  return value;
}

float Utilities::ReturnPhi(const reco::Candidate::LorentzVector& theA, const reco::Candidate::LorentzVector& theL1, const reco::Candidate::LorentzVector& theL2, const reco::Candidate::LorentzVector& theB1, const reco::Candidate::LorentzVector& theB2) {
  TLorentzVector pA(theA.px(), theA.py(), theA.pz(), theA.energy());
  TLorentzVector pL1(theL1.px(), theL1.py(), theL1.pz(), theL1.energy());
  TLorentzVector pL2(theL2.px(), theL2.py(), theL2.pz(), theL2.energy());
  TLorentzVector pB1(theB1.px(), theB1.py(), theB1.pz(), theB1.energy());
  TLorentzVector pB2(theB2.px(), theB2.py(), theB2.pz(), theB2.energy());
  // Boost objects to the A rest frame
  pL1.Boost( -pA.BoostVector() );
  pL2.Boost( -pA.BoostVector() );
  pB1.Boost( -pA.BoostVector() );
  pB2.Boost( -pA.BoostVector() );
  // Build unit vectors orthogonal to the decay planes
  TVector3 Zplane=pL1.Vect().Cross( pL2.Vect() ); // L1 x L2
  TVector3 Hplane=pB1.Vect().Cross( pB2.Vect() ); // B1 x B2
  Zplane.SetMag(1.);
  Hplane.SetMag(1.);
  // Sign of Phi
  TLorentzVector pZr = pL1 + pL2;
  float sgn = pZr.Vect().Dot( Zplane.Cross(Hplane) );
  sgn/=fabs(sgn);
  
  float value=sgn * acos( Zplane.Dot(Hplane) );
  if(value!=value || isinf(value)) return 0.;
  return value;
}


float Utilities::ReturnPhi1(const reco::Candidate::LorentzVector& theA, const reco::Candidate::LorentzVector& theL1, const reco::Candidate::LorentzVector& theL2) {
  TLorentzVector pA(theA.px(), theA.py(), theA.pz(), theA.energy());
  TLorentzVector pL1(theL1.px(), theL1.py(), theL1.pz(), theL1.energy());
  TLorentzVector pL2(theL2.px(), theL2.py(), theL2.pz(), theL2.energy());
  TVector3 beamAxis(0., 0., 1.);
  // Boost objects to the A rest frame
  pL1.Boost( -pA.BoostVector() );
  pL2.Boost( -pA.BoostVector() );
  // Reconstruct Z in A rest frame
  TLorentzVector pZr = pL1 + pL2;
  // Build unit vectors orthogonal to the decay planes
  TVector3 Zplane=pL1.Vect().Cross( pL2.Vect() ); // L1 x L2
  TVector3 Bplane=beamAxis.Cross( pZr.Vect() ); // Beam x Z, beam/Z plane
  Zplane.SetMag(1.);
  Bplane.SetMag(1.);
  // Sign of Phi1
  float sgn = pZr.Vect().Dot( Zplane.Cross(Bplane) );
  sgn/=fabs(sgn);
  
  float value=sgn * acos( Zplane.Dot(Bplane) );
  if(value!=value || isinf(value)) return 0.;
  return value;
}


// Kinematics and reconstruction


// -----------------------------------
// ---------- KINEMATIC FIT ----------
// -----------------------------------

float Utilities::PerformKinematicFit(pat::Jet* tJet1, pat::Jet* tJet2, reco::Candidate::LorentzVector* fJet1, reco::Candidate::LorentzVector* fJet2, float mass) {
  //  TLorentzVector b1, b2;
  //  b1.SetPtEtaPhiE(tJet1->pt(), tJet1->eta(), tJet1->phi(), tJet1->energy());
  //  b2.SetPtEtaPhiE(tJet2->pt(), tJet2->eta(), tJet2->phi(), tJet2->energy());
    
    TMatrixD m1(3,3);
    TMatrixD m2(3,3);
    m1.Zero();
    m2.Zero();

    //In this example the covariant matrix depends on the transverse energy and eta of the jets
    m1(0,0) = GetErrEt (tJet1->et(), tJet1->eta()); // et
    m1(1,1) = GetErrEta(tJet1->et(), tJet1->eta()); // eta
    m1(2,2) = GetErrPhi(tJet1->et(), tJet1->eta()); // phi
    m2(0,0) = GetErrEt (tJet2->et(), tJet2->eta()); // et
    m2(1,1) = GetErrEta(tJet2->et(), tJet2->eta()); // eta
    m2(2,2) = GetErrPhi(tJet2->et(), tJet2->eta()); // phi

  //  TFitParticleEtEtaPhi jet1("Jet1", "Jet1", &b1, &m1);
  //  TFitParticleEtEtaPhi jet2("Jet2", "Jet2", &b2, &m2);
    
  //  TVector3 b1_3=b1.Vect();
  //  TVector3 b2_3=b2.Vect();
    TVector3 b1(tJet1->px(), tJet1->py(), tJet1->pz());
    TVector3 b2(tJet2->px(), tJet2->py(), tJet2->pz());
    TFitParticlePtEtaPhi jet1("Jet1", "Jet1", &b1, tJet1->mass(), &m1 );
    TFitParticlePtEtaPhi jet2("Jet2", "Jet2", &b2, tJet2->mass(), &m2 );

  //  TFitParticleEScaledMomDev jet1("Jet1", "Jet1", &b1, &m1);
  //  TFitParticleEScaledMomDev jet2("Jet2", "Jet2", &b2, &m2);

    //vec1 and vec2 must make a W boson
    TFitConstraintM mCons1("hMassConstraint", "hMass-Constraint", 0, 0, mass);
    mCons1.addParticles1( &jet1, &jet2 );

    //Definition of the fitter
    //Add two constraints
    TKinFitter fitter("fitter", "fitter");
    fitter.addMeasParticle( &jet1 );
    fitter.addMeasParticle( &jet2 );

    fitter.addConstraint( &mCons1 );
    
    //Set convergence criteria
    fitter.setMaxNbIter( 30 );
    fitter.setMaxDeltaS( 1e-2 );
    fitter.setMaxF( 1e-1 );
    fitter.setVerbosity(1);

    // Perform the fit
    fitter.fit();
  //  fitter.print();
    
    
    float dPt1  = jet1.getCurr4Vec()->Pt()  - jet1.getIni4Vec()->Pt();
    float dEta1 = jet1.getCurr4Vec()->Eta() - jet1.getIni4Vec()->Eta();
    float dPhi1 = jet1.getCurr4Vec()->Phi() - jet1.getIni4Vec()->Phi();
    float dPt2  = jet2.getCurr4Vec()->Pt()  - jet2.getIni4Vec()->Pt();
    float dEta2 = jet2.getCurr4Vec()->Eta() - jet2.getIni4Vec()->Eta();
    float dPhi2 = jet2.getCurr4Vec()->Phi() - jet2.getIni4Vec()->Phi();

    float chi2( dPt1*dPt1/m1(0,0) + dEta1*dEta1/m1(1,1) + dPhi1*dPhi1/m1(2,2) + dPt2*dPt2/m2(0,0) + dEta2*dEta2/m2(1,1) + dPhi2*dPhi2/m2(2,2) ); //=fitter.getS();
    /*
    float pchi2=TMath::Prob(chi2, fitter.getNDF());
    
    Hist["k_chi2"]->Fill(chi2, EventWeight);
    Hist["k_chi2Prob"]->Fill(pchi2, EventWeight);
    Hist["k_deltaPt1" ]->Fill(dPt1, EventWeight);
    Hist["k_deltaEta1"]->Fill(dEta1, EventWeight);
    Hist["k_deltaPhi1"]->Fill(dPhi1, EventWeight);
    Hist["k_deltaPt2" ]->Fill(dPt2, EventWeight);
    Hist["k_deltaEta2"]->Fill(dEta2, EventWeight);
    Hist["k_deltaPhi2"]->Fill(dPhi2, EventWeight);
    if(tJet1->genParton()) {
      Hist["k_pullPt1" ]->Fill((jet1.getCurr4Vec()->Pt()-tJet1->genParton()->pt())/sqrt(m1(0,0)), EventWeight);
      Hist["k_pullEta1"]->Fill((jet1.getCurr4Vec()->Eta()-tJet1->genParton()->eta())/sqrt(m1(1,1)), EventWeight);
      Hist["k_pullPhi1"]->Fill((jet1.getCurr4Vec()->Phi()-tJet1->genParton()->phi())/sqrt(m1(2,2)), EventWeight);
    }
    if(tJet2->genParton()) {
      Hist["k_pullPt2" ]->Fill((jet2.getCurr4Vec()->Pt()-tJet2->genParton()->pt())/sqrt(m2(0,0)), EventWeight);
      Hist["k_pullEta2"]->Fill((jet2.getCurr4Vec()->Eta()-tJet2->genParton()->eta())/sqrt(m2(1,1)), EventWeight);
      Hist["k_pullPhi2"]->Fill((jet2.getCurr4Vec()->Phi()-tJet2->genParton()->phi())/sqrt(m2(2,2)), EventWeight);
    }
    */
    // Update objects
    fJet1->SetPxPyPzE(jet1.getCurr4Vec()->Px(), jet1.getCurr4Vec()->Py(), jet1.getCurr4Vec()->Pz(), jet1.getCurr4Vec()->Energy());
    fJet2->SetPxPyPzE(jet2.getCurr4Vec()->Px(), jet2.getCurr4Vec()->Py(), jet2.getCurr4Vec()->Pz(), jet2.getCurr4Vec()->Energy());
    
    return chi2;
}



// Neutrino pz recovering by trading W mass

float Utilities::RecoverNeutrinoPz(const reco::Particle::LorentzVector* lep, const reco::Particle::LorentzVector* met) {
    // W kinematical reconstruction
    float pz = 0.;
    float a = pow(80.4,2) - pow(lep->mass(),2) + 2.*lep->px()*met->px() + 2.*lep->py()*met->py();
    float A = 4*( pow(lep->energy(),2) - pow(lep->pz(),2) );
    float B = -4*a*lep->pz();
    float C = 4*pow(lep->energy(),2) * (pow(met->px(),2)  + pow(met->py(),2)) - pow(a,2);
    float D = pow(B,2) - 4*A*C;
    // If there are real solutions, use the one with lowest pz                                            
    if (D>=0) {
        float s1 = (-B+sqrt(D))/(2*A);
        float s2 = (-B-sqrt(D))/(2*A);
        if(fabs(s1)<fabs(s2)) pz=s1;
        else pz=s2;
    }
    // Otherwise, use real part                                                                           
    else {
        pz = -B/(2*A);
    }
    return pz;
}






pat::CompositeCandidate Utilities::RecoilMassFormula(pat::CompositeCandidate& H, pat::MET& met){
    pat::CompositeCandidate X;
    AddFourMomenta addP4;
    X.addDaughter(H);
    X.addDaughter(met);
    addP4.set(X);
    reco::Particle::LorentzVector metp4 = met.p4();
    reco::Particle::LorentzVector Xp4;
    metp4.SetPz(-H.pz());
    Xp4 += metp4;
    Xp4.SetPz(0);
    float B = -2.*H.energy();
    float C = pow(H.mass(),2) - pow(90.18,2);
    float D = pow(B,2) - 4*1*C;
    float mX;
    if(D>0){
        float s1 = (-B+sqrt(D))/2.;
        float s2 = (-B-sqrt(D))/2.;
        if(fabs(s1)>fabs(s2)) mX = s1;
        else mX = s2;
    }
    else{
        mX = -B/2.;
    }
    Xp4.SetE(sqrt(pow(mX,2) + pow(Xp4.Px(),2) + pow(Xp4.Py(),2) + pow(Xp4.Pz(),2)));
    X.setP4(Xp4);
    X.setCharge(0);
    return X;
}

//Find mother id of a const reco::GenParticle, method used for genPartons MC truth
int Utilities::FindMotherId(const reco::GenParticle* p) {
  int pId = p->pdgId();
  const reco::Candidate* mom = p->mother();
  while (mom != 0 && mom->pdgId() == pId)
    mom = mom->mother();
  return mom->pdgId();
}

#endif
