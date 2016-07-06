#ifndef PROT_FEATURE_PREC_ENV_HPP_
#define PROT_FEATURE_PREC_ENV_HPP_

#include <memory>
#include <vector>

#include "feature/feature_mng.hpp" 
#include "feature/real_env.hpp" 

namespace prot {

class PrecEnv {
 public:
  static RealEnvPtr deconv(double prec_win_size, PeakPtrVec &peak_list, 
                           double prec_mz, int prec_charge);
};


}

#endif
