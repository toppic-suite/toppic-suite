#ifndef PROT_FEATURE_DECONV_ONE_SP_HPP_
#define PROT_FEATURE_DECONV_ONE_SP_HPP_

#include "feature/feature_mng.hpp"
#include "feature/deconv_data.hpp"
#include "feature/match_env.hpp"

namespace prot {

class DeconvOneSp {
 public:
  DeconvOneSp(FeatureMngPtr mng_ptr);

  void setData(PeakPtrVec &peak_list);

  void preprocess();

 private:
	FeatureMngPtr mng_ptr_;
	DeconvDataPtr data_ptr_;
	MatchEnvPtrVec result_envs_;
};


}
#endif
