#ifndef PROT_FEATURE_ENV_DETECT_HPP_
#define PROT_FEATURE_ENV_DETECT_HPP_

#include <memory>
#include <vector>

#include "feature/envelope.hpp"
#include "feature/match_env.hpp"
#include "feature/feature_mng.hpp"
#include "feature/deconv_data.hpp"

namespace prot {

class EnvDetect {
 public:
  static double calcInteRatio(EnvelopePtr theo_env, PeakPtrVec &peak_list, 
                              double tolerance);
  static MatchEnvPtr detectEnv(PeakPtrVec &peak_list, int base_peak,
                               int charge, double max_mass, FeatureMngPtr mng_ptr);

  static MatchEnvPtr detectEnv(PeakPtrVec &peak_list, double mono_mass, 
                               int charge, FeatureMngPtr mng_ptr);

  static MatchEnvPtr2D getCandidate(DeconvDataPtr data_ptr, FeatureMngPtr mng_ptr);
};

}

#endif
