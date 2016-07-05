#ifndef PROT_FEATURE_CHARGE_CMP_HPP_
#define PROT_FEATURE_CHARGE_CMP_HPP_

#include <memory>
#include <vector>

#include "spec/peak.hpp"
#include "feature/match_env.hpp"

namespace prot {

class ChargeCmp {
 public:
  static int comp(PeakPtrVec &peak_list, MatchEnvPtr a, MatchEnvPtr b, double tolerance);
};

}

#endif
