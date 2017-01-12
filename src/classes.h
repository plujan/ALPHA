
#include <vector>

#include "Analysis/ALPHA/plugins/Objects.h"
#include "Analysis/ALPHA/interface/alp_objects.h"

namespace { 
  struct dictionary {
    CandidateType dummy0;
    JetType dummy1;
    LeptonType dummy2;
    std::vector<CandidateType> dummy3;
    std::vector<JetType> dummy4;
    std::vector<LeptonType> dummy5;
    std::vector<std::size_t> dummy6;
    LorentzType dummy7;
    std::vector<LorentzType> dummy8;

    alp::Candidate alp_candidate_;
    std::vector<alp::Candidate> vector_alp_candidate_;
    alp::Jet alp_jet_;
    std::vector<alp::Jet> vector_alp_jet_;
    alp::DiObject alp_diobject_;
    std::vector<alp::DiObject> vector_alp_diobject_;


  };
}
