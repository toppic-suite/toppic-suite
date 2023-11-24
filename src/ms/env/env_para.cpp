//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include "common/util/logger.hpp"
#include "ms/env/env_para.hpp" 

namespace toppic {

EnvPara::EnvPara(double mz_tolerance) {
  mz_tolerance_ = mz_tolerance; 
  score_error_tolerance_ = mz_tolerance; 
}

// get the mass group based on mass value 
int EnvPara::getMassGroup(double base_mass) {
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
int EnvPara::compMinConsPeakNum(int peak_num, int mass_group) {
  int min_cons_peak_num = peak_num + relative_consecutive_peak_num_;
  if (min_cons_peak_num < min_consecutive_peak_num_[mass_group]) {
    min_cons_peak_num = min_consecutive_peak_num_[mass_group];
  }
  return min_cons_peak_num;
}

int EnvPara::compLowMassNum() {
  int low_mass_num = static_cast<int>(low_high_dividor_ / aa_avg_mass_ * peak_density_);
  return low_mass_num;
}

int EnvPara::compHighMassNum(double prec_mass) {
  int high_mass_num = static_cast<int> ((prec_mass - low_high_dividor_) / aa_avg_mass_ * peak_density_);
  if (high_mass_num < 0) {
    high_mass_num = 0;
  }
  return high_mass_num;
}

}
