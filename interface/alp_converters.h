
#pragma once

#include <vector>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"

#include "alp_objects.h"

namespace alp {

  void convert(const pat::Jet & orig, alp::Jet & dest) {

    dest.p4_ = orig.p4();

    // generator/MC inforamtion
    if (orig.genJet()) {
      dest.genJets_.emplace_back(orig.genJet()->p4());
    }
    if (orig.genParton()) {
      dest.genPartons_.emplace_back(orig.genParton()->p4());
    }
    dest.partonFlavour_ = orig.partonFlavour();
    dest.hadronFlavour_ = orig.hadronFlavour();
    dest.ptRaw_ = orig.correctedJet(0).pt();
    dest.JESunc_ = orig.userFloat("JESUncertainty");

    // bDiscriminators
    std::vector<std::string> bDiscs = { "pfCombinedInclusiveSecondaryVertexV2BJetTags",
                                        "pfCombinedMVAV2BJetTags" };
    for (const auto & bDisc : bDiscs) {
      dest.discs_.emplace_back(bDisc, orig.bDiscriminator(bDisc));
    }
    // other discriminators (userFloats) 
    std::vector<std::string> oDiscs = { "ReshapedDiscriminator",
                                        "ReshapedDiscriminatorUp",
                                        "ReshapedDiscriminatorDown",
                                        "CMVAR",
                                        "CMVARUp",
                                        "CMVARDown",
                                        "QGLikelihood"
                                       };
    for (const auto & oDisc : oDiscs) {
      if (orig.hasUserFloat(oDisc)) {
        dest.discs_.emplace_back(oDisc, orig.userFloat(oDisc));
      }
    }
  }


  void convert(const std::vector<pat::Jet> & orig, std::vector<alp::Jet> & dest) {
    dest.clear();
    for (const auto & e : orig) {
      dest.emplace_back();
      convert(e,dest.back());
    }
  }

}
  

