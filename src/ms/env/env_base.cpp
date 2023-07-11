//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <fstream>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "ms/env/env_base.hpp"

namespace toppic {

EnvBasePtr EnvBase::env_base_ptr_;

EnvBase::EnvBase(std::string file_name, int entry_num, 
                 double mass_interval):
  entry_num_(entry_num),
  mass_interval_(mass_interval) {
    std::ifstream input;
    input.open(file_name.c_str(), std::ios::in);
    if (!input.is_open()) {
      LOG_ERROR("Env file  " << file_name << " does not exist.");
      exit(EXIT_FAILURE); 
    }
    for (int i = 0; i < entry_num_; i++) {
      int peak_num = 0;
      std::string line;
      std::vector<std::string> line_list;
      while (std::getline(input, line)) {
        if (line == "") {
          break;
        }
        peak_num++;
        line_list.push_back(line);
      }
      envs_.push_back(std::make_shared<Env>(peak_num - 1, line_list));
    }
    input.close();
    initRefMassIdx();
  }

void EnvBase::initBase(const std::string &resource_dir) {
  std::string distr_file_name = resource_dir + file_util::getFileSeparator() + 
      getBaseDirName() + file_util::getFileSeparator() + getBaseFileName();
  LOG_DEBUG("distribution file name " << distr_file_name);
  int distr_entry_num = getDistrEntryNum();
  double distr_mass_interval = getDistrMassInterval();
  env_base_ptr_ = std::make_shared<EnvBase>(distr_file_name, distr_entry_num, distr_mass_interval);
}

void EnvBase::initRefMassIdx() {
  ref_mass_idxes_.resize(entry_num_, -1);
  for (int i = 0; i < entry_num_; i++) {
    double mz = envs_[i]->getReferMz();
    int idx = static_cast<int>(mz / mass_interval_);
    if (idx >= 0 && idx < entry_num_) {
     ref_mass_idxes_[idx] = i;
    }
  }
  ref_mass_idxes_[0] = 0;
  for (int i = 1; i < entry_num_; i++) {
    if (ref_mass_idxes_[i] < 0) {
      ref_mass_idxes_[i] = ref_mass_idxes_[i - 1];
    }
  }
}

// Private class methods
// Get a distribution envelope based on its monoisotopic mass.
EnvPtr EnvBase::getBaseEnvByMonoMass(double mass) {
  int idx = static_cast<int>(mass / mass_interval_);
  if (idx < 0) {
    LOG_ERROR("Invalid mass!");
    exit(EXIT_FAILURE);
  } else if (idx >= entry_num_) {
    LOG_WARN("Mass out of bound!");
    return nullptr;
  }
  return envs_[idx];
}

// Get a distribution envelope based on its mass of reference (highest
// intensity) peak.
EnvPtr EnvBase::getBaseEnvByRefMass(double mass) {
  int idx = static_cast<int>(mass / mass_interval_);
  if (idx < 0) {
    LOG_ERROR("Invalid mass");
    exit(EXIT_FAILURE);
  } else if (idx >= entry_num_) {
    LOG_WARN("Mass out of bound!");
    return nullptr;
  }
  return envs_[ref_mass_idxes_[idx]];
}


// public static methods
EnvPtr EnvBase::getEnvByMonoMass(double mass) {
  return env_base_ptr_->getBaseEnvByMonoMass(mass);
}

EnvPtr EnvBase::getEnvByRefMass(double mass) {
  return env_base_ptr_->getBaseEnvByRefMass(mass);
}

double EnvBase::convertMonoMassToAvgMass(double mass) {
  EnvPtr env_ptr = env_base_ptr_->getBaseEnvByMonoMass(mass);
  if (env_ptr == nullptr) {
    LOG_ERROR("Invalid mass!");
    exit(EXIT_FAILURE);
  }
  double diff = env_ptr->getAvgNeutralMass() - env_ptr->getMonoNeutralMass();
  return mass + diff;
}

//Ref mass is the monoisotopic mass of the highest peak in the envelope
double EnvBase::convertRefMassToMonoMass(double mass) {
  EnvPtr env_ptr = env_base_ptr_->getBaseEnvByRefMass(mass);
  if (env_ptr == nullptr) {
    LOG_ERROR("Invalid mass!");
    exit(EXIT_FAILURE);
  }
  double diff = env_ptr->getReferNeutralMass() - env_ptr->getMonoNeutralMass();
  double mono_mass = mass - diff;
  if (mono_mass < 0) {
    mono_mass = 0;
  }
  return mono_mass; 
}

}  // namespace toppic
