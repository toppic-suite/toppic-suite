//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include <string>
#include <vector>

#include "base/logger.hpp"
#include "feature/env_base.hpp"

namespace prot {

EnvBase::EnvBase(std::string file_name, int entry_num, double mass_interval):
    entry_num_(entry_num),
    mass_interval_(mass_interval) {
      std::ifstream input;
      input.open(file_name.c_str(), std::ios::in);
      if (!input.is_open()) {
        LOG_ERROR("env file  " << file_name << " does not exist.");
        throw "env file does not exist.";
      }
      LOG_DEBUG("start reading");
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
        // LOG_DEBUG("one envelope");
        envs_.push_back(std::make_shared<Envelope>(peak_num - 1, line_list));
      }
      input.close();
      initBaseMassIdx();
    }

void EnvBase::initBaseMassIdx() {
  base_mass_idxes_.resize(entry_num_, -1);
  for (int i = 0; i < entry_num_; i++) {
    double mz = envs_[i]->getReferMz();
    int idx = static_cast<int>(mz / mass_interval_);
    if (idx >= 0 && idx < entry_num_) {
      base_mass_idxes_[idx] = i;
    }
  }
  base_mass_idxes_[0] = 0;
  for (int i = 1; i < entry_num_; i++) {
    if (base_mass_idxes_[i] < 0) {
      base_mass_idxes_[i] = base_mass_idxes_[i - 1];
    }
  }
}

// Get a distribution envelope based on its monoisotopic mass.
EnvelopePtr EnvBase::getEnvByMonoMass(double mass) {
  int idx = static_cast<int>(mass / mass_interval_);
  if (idx < 0) {
    LOG_ERROR("Invalid mass");
    exit(1);
  } else if (idx >= entry_num_) {
    LOG_DEBUG("mass out of bound");
    return nullptr;
  }
  return envs_[idx];
}

// Get a distribution envelope based on its mass of base peak.
EnvelopePtr EnvBase::getEnvByBaseMass(double mass) {
  int idx = static_cast<int>(mass / mass_interval_);
  if (idx < 0) {
    LOG_ERROR("Invalid mass");
    exit(1);
  } else if (idx >= entry_num_) {
    LOG_DEBUG("Mass out of bound");
    return nullptr;
  }
  return envs_[base_mass_idxes_[idx]];
}

}  // namespace prot
