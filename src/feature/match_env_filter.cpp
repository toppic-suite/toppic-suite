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
