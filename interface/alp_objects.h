
#pragma once

#include <vector>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"

namespace alp {
  
  typedef ROOT::Math::PtEtaPhiEVector PtEtaPhiEVector;
  typedef std::pair<std::string, float> StringFloatPair;
  typedef std::vector<StringFloatPair> StringFloatPairVector;
  typedef std::pair<std::string, int> StringIntPair;
  typedef std::vector<StringIntPair> StringIntPairVector;

  constexpr auto CSV_name = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
  constexpr auto CMVA_name = "pfCombinedMVAV2BJetTags";


  class Candidate {

    public:

      // default constructor
      Candidate() :
        p4_() {}
      // copy constructor
      Candidate( const Candidate& rhs) : 
        p4_(rhs.p4_) {}
      // constructors from four vector
      Candidate( const PtEtaPhiEVector & p4) :
        p4_(p4) {}

      double pt() const {return p4_.Pt();}
      double eta() const {return p4_.Eta();}
      double phi() const {return p4_.Phi();}
      double energy() const {return p4_.E();}
      double et() const {return p4_.Et();}

    // attributes  (also public) 
      PtEtaPhiEVector p4_;
  };

  typedef std::vector<alp::Candidate> CandidateCollection;

  class Jet : public Candidate {

    public:
      // default constructor
      Jet() : Candidate() {}
      // copy constructor
      Jet( const Jet & rhs) : 
        Candidate(rhs),
        partonFlavour_(rhs.partonFlavour_),
        hadronFlavour_(rhs.hadronFlavour_),
        ptRaw_(rhs.ptRaw_),
        JESunc_(rhs.JESunc_),
        discs_(rhs.discs_),
        ids_(rhs.ids_) {}

      // inherit other constructors
      using Candidate::Candidate;

      int id(const std::string & name ) const {
        for( auto id : ids_)  {
          if (id.first == name) return id.second;
        }
        return -777; 
      }

      int idC(const char * name ) const {
        return id(std::string(name));
      }

      float disc(const std::string & name ) const {
        for( auto disc : discs_)  {
          if (disc.first == name) return disc.second;
        }
        return -777.; 
      }

      float discC(const char * name ) const {
        return disc(std::string(name));
      }

      int partonFlavour() const {return partonFlavour_;};
      int hadronFlavour() const {return hadronFlavour_;};
      float CSV() const { return disc(CSV_name);}
      float CMVA() const { return disc(CMVA_name);}
      float ptRaw() const { return ptRaw_;};
      float JESunc() const { return JESunc_;};
      
    // attributes (also public)

      // flavour attributes (0 if data/undefined)
      int partonFlavour_ = 0;
      int hadronFlavour_ = 0;
      float ptRaw_ = 0;
      float JESunc_ = 0;

      // to keep float discriminator values
      StringFloatPairVector discs_ = {};
      // to keep integer ids
      StringIntPairVector ids_ = {};

      std::vector<PtEtaPhiEVector> genJets_;
      std::vector<PtEtaPhiEVector> genPartons_;

  };

  typedef std::vector<alp::Jet> JetCollection;

}
