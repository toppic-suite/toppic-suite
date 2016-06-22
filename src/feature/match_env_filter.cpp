#include <algorithm>

#include "base/logger.hpp"
#include "feature/match_env_filter.hpp" 

namespace prot {

MatchEnvPtrVec MatchEnvFilter::filter(MatchEnvPtrVec &ori_envs, double prec_mass,
                                      FeatureMngPtr mng_ptr) {
  MatchEnvPtrVec low_mass_envs;
  MatchEnvPtrVec high_mass_envs;
  std::sort(ori_envs.begin(), ori_envs.end(), MatchEnv::cmpScoreDec);
  int low_mass_num = (int) (mng_ptr->low_high_dividor_ / mng_ptr->aa_avg_mass_ * mng_ptr->peak_density_);
  int high_mass_num = (int) ((prec_mass - mng_ptr->low_high_dividor_)
                         / mng_ptr->aa_avg_mass_ * mng_ptr->peak_density_);
  for (size_t i = 0; i < ori_envs.size(); i++) {
    if (ori_envs[i]->getRealEnvPtr()->getMonoMass() <= mng_ptr->low_high_dividor_) {
      if ((int)low_mass_envs.size() < low_mass_num) {
        low_mass_envs.push_back(ori_envs[i]);
      }
    } else {
      if ((int)high_mass_envs.size() < high_mass_num) {
        high_mass_envs.push_back(ori_envs[i]);
      }
    }
  }
  MatchEnvPtrVec result;
  result.insert(std::end(result), std::begin(low_mass_envs), std::end(low_mass_envs));
  result.insert(std::end(result), std::begin(high_mass_envs), std::end(high_mass_envs));
  std::sort(result.begin(), result.end(), MatchEnv::cmpScoreDec); 
  return result;
}

}
