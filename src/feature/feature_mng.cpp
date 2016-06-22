#include "feature/feature_mng.hpp" 

namespace prot {

FeatureMng::FeatureMng() {
  env_base_ptr_ = EnvBasePtr(new EnvBase(distr_file_name_, distr_entry_num_, distr_mass_interval_));
}

void FeatureMng::setMinInte(double min_inte) {
  min_inte_ = min_inte;
  min_refer_inte_ = min_inte * sn_ratio_;
}

// get the mass group based on mass value 
int FeatureMng::getMassGroup(double base_mass) {
  int group = -1;
  for (int i = 0; i < (int)mass_group_boundary_.size() - 1; i++) {
    if (base_mass >= mass_group_boundary_[i] && base_mass < mass_group_boundary_[i + 1]) {
      group = i;
      break;
    }
  }
  return group;
}

// compute minimum consecutive peak num: max{peak_num -3 , 3} 
int FeatureMng::compMinConsPeakNum(int peak_num, int mass_group) {
  int min_cons_peak_num = peak_num + relative_consecutive_peak_num_;
  if (min_cons_peak_num < min_consecutive_peak_num_[mass_group]) {
    min_cons_peak_num = min_consecutive_peak_num_[mass_group];
  }
  return min_cons_peak_num;
}

void FeatureMng::setTolerance(double tolerance) {
  mz_tolerance_ = tolerance;
  score_error_tolerance_ = tolerance;
}

}
