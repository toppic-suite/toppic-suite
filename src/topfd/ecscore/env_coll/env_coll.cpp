//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
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
#include "topfd/ecscore/env_coll/env_coll.hpp"

namespace toppic {

EnvColl::EnvColl(SeedEnvelopePtr seed_ptr, EnvSetPtrVec &env_set_list,
                 int min_charge, int max_charge, 
                 int start_spec_id, int end_spec_id) {
  seed_ptr_ = seed_ptr;
  env_set_list_ = env_set_list;
  min_charge_ = min_charge;
  max_charge_ = max_charge;
  start_spec_id_ = start_spec_id;
  end_spec_id_ = end_spec_id;
  exp_inte_sum_list_ = compExpInteSumList();
}

std::vector<double> EnvColl::compExpInteSumList() {
  int peak_num = seed_ptr_->getPeakNum();
  std::vector<double> sum_list(peak_num, 0);
  for (auto &env_set: env_set_list_) {
    if (env_set == nullptr) {
      continue;
    }
    std::vector<double> cur_sum_list = env_set->compExpInteSumList(); 
    if (peak_num != static_cast<int>(cur_sum_list.size())) {
      LOG_ERROR("peak number " << peak_num << " cur list len " << cur_sum_list.size());
    }
    for (int i = 0; i < peak_num; i++)
      sum_list[i] = sum_list[i] + cur_sum_list[i];
  }
  return sum_list;
}

std::vector<int> EnvColl::getChargeList() {
  std::vector<int> charge_list;
  for (auto es: env_set_list_)
    charge_list.push_back(es->getCharge());
  return charge_list;
}

void EnvColl::refineMonoMass() {
  double weight = 0;
  double weight_mz_error = 0;
  for (auto &env_set: env_set_list_) {
    double cur_weight = 0, cur_weight_mz_error = 0;
    env_set->getWeightMzError(cur_weight, cur_weight_mz_error);
    weight = weight + cur_weight;
    weight_mz_error = weight_mz_error + cur_weight_mz_error;
  }
  if (weight > 0) {
    double mz_error = weight_mz_error / weight;
    seed_ptr_->shift(mz_error * seed_ptr_->getCharge());
  }
  else {
    LOG_ERROR("ERROR 0 weight in refine_mono_mass");
  }
}

double EnvColl::getIntensity(double sn_ratio, double noise_inte) {
  double inte = 0;
  for (auto env_set: env_set_list_) {
    double tmp_inte = env_set->compIntensity(sn_ratio, noise_inte);
    inte = inte + tmp_inte;
  }
  return inte;
}

EnvSetPtr EnvColl::getSeedEnvSet() {
  for (const auto &es: env_set_list_) {
    SeedEnvelopePtr es_seed_env = es->getSeedPtr();
    if (es_seed_env->getCharge() == seed_ptr_->getCharge()) {
      return es;
    }
  }
  return nullptr;
}

}
