// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "base/logger.hpp"
#include "feature/env_base.hpp" 

namespace prot {

EnvBase::EnvBase(std::string file_name, int entry_num, 
                 double mass_interval):
    entry_num_(entry_num),
    mass_interval_(mass_interval) {
      std::ifstream input;
      input.open(file_name.c_str(), std::ios::in);
      if (!input.is_open()) {
        LOG_ERROR( "env file  " << file_name << " does not exist.");
        throw "fasta file does not exist.";
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
        //LOG_DEBUG("one envelope");
        envs_.push_back(EnvelopePtr(new Envelope(peak_num - 1, line_list)));
      }
      input.close();
      initBaseMassIdx();
    }

void EnvBase::initBaseMassIdx() {
  base_mass_idxes_.resize(entry_num_, -1);
  for (int i = 0; i < entry_num_; i++) {
    double mz = envs_[i]->getReferMz();
    int idx = (int) (mz / mass_interval_);
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
  int idx = (int) (mass / mass_interval_);
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
  int idx = (int) (mass / mass_interval_);
  if (idx < 0) {
    LOG_ERROR("Invalid mass");
    exit(1);
  } else if (idx >= entry_num_) {
    LOG_DEBUG("Mass out of bound");
    return nullptr;
  }
  return envs_[base_mass_idxes_[idx]];
}

}
