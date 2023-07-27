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

#include <algorithm>

#include "common/util/logger.hpp"
#include "topfd/ecscore/env_coll/env_coll.hpp"

namespace toppic {

EnvColl::EnvColl(SeedEnvPtr seed_ptr, EnvSetPtrVec &env_set_list,
                 int min_charge, int max_charge,
                 int start_spec_id, int end_spec_id) {
  seed_ptr_ = seed_ptr;
  env_set_list_ = env_set_list;
  min_charge_ = min_charge;
  max_charge_ = max_charge;
  start_spec_id_ = start_spec_id;
  end_spec_id_ = end_spec_id;
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
    std::pair<double, double> error_weight = env_set->getMzErrorAndWeight();
    weight_mz_error = weight_mz_error + error_weight.first;
    weight = weight + error_weight.second;
  }
  if (weight > 0) {
    double mz_error = weight_mz_error / weight;
    seed_ptr_->changeMz(mz_error); 
  }
  else {
    LOG_INFO("ERROR 0 weight in refine_mono_mass");
  }
}

double EnvColl::getIntensity() {
  double inte = 0;
  for (auto env_set: env_set_list_) {
    double tmp_inte = env_set->getInte();
    inte = inte + tmp_inte;
  }
  return inte;
}

EnvSetPtr EnvColl::getSeedEnvSet() {
  for (const auto &es: env_set_list_) {
    SeedEnvPtr es_seed_env = es->getSeedPtr();
    if (es_seed_env->getCharge() == seed_ptr_->getCharge()) {
      return es;
    }
  }
  return nullptr;
}

void EnvColl::removePeakData(MsMapPtr matrix_ptr) {
  for (auto env_set_ptr: env_set_list_) {
    if (env_set_ptr != nullptr) {
      env_set_ptr->removePeakData(matrix_ptr); 
    }
  }
}

void EnvColl::mergeEnvSet(EnvSetPtr new_set_ptr) {
  new_set_ptr->setSeedPtr(seed_ptr_);
  if (new_set_ptr->getStartSpecId() < start_spec_id_) {
    start_spec_id_ = new_set_ptr->getStartSpecId();
  }
  if (new_set_ptr->getEndSpecId() > end_spec_id_) {
    end_spec_id_ = new_set_ptr->getEndSpecId();
  }
  int charge = new_set_ptr->getCharge();
  for (size_t i = 0; i < env_set_list_.size(); i++) {
    if (env_set_list_[i]->getCharge() == charge) {
      env_set_list_[i]->mergeEnvSet(new_set_ptr);
      return;
    }
  }
  // if no matched charge is found, add it to the env_set_list;
  if (charge < min_charge_) {
    min_charge_ = charge;
  }
  if (charge > max_charge_) {
    max_charge_ = charge;
  }
  env_set_list_.push_back(new_set_ptr);
  std::sort(env_set_list_.begin(), env_set_list_.end(),EnvSet::cmpChargeInc); 
}

}
