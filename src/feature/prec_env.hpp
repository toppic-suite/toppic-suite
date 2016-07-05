#ifndef PROT_FEATURE_PREC_ENV_HPP_
#define PROT_FEATURE_PREC_ENV_HPP_

#include <memory>
#include <vector>

#include "feature/feature_mng.hpp" 
#include "feature/real_env.hpp" 

namespace prot {

struct PrecInfo {
  double mono_mz = 0;
  double avg_mz = 0;
  int charge = 0;
};

class PrecEnv {
 public:
  static RealEnvPtr deconv(double prec_win_size, PeakPtrVec &peak_list, 
                           double prec_mz, double prec_charge);
};


}

#endif
