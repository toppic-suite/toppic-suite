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

#ifndef TOPPIC_TOPFD_ECSCORE_ENV_SEED_ENV_HPP
#define TOPPIC_TOPFD_ECSCORE_ENV_SEED_ENV_HPP

#include <vector>

#include "ms/spec/deconv_ms.hpp"
#include "ms/env/env_base.hpp"

namespace toppic {

class SeedEnv;

typedef std::shared_ptr<SeedEnv> SeedEnvPtr;

class SeedEnv : public Env {
 public:
  SeedEnv(DeconvPeakPtr peak_ptr);

  SeedEnv(SeedEnvPtr env_ptr, int new_charge);

  int getSpecId() const { return spec_id_; }

  void setSpecId(int spec_id) { spec_id_ = spec_id; }

  double getSeedInte() const {return seed_inte_;}

  static bool cmpSeedInteDec(const SeedEnvPtr a, const SeedEnvPtr b) {
    return a->getSeedInte() > b->getSeedInte(); }

  std::string getString();

  bool containEnoughPeaks();

 private:
  int spec_id_;
  double seed_inte_;
};

typedef std::vector<SeedEnvPtr> SeedEnvPtrVec;
typedef std::vector<SeedEnvPtrVec> SeedEnvPtr2D;

}

#endif //TOPPIC_SEED_ENVELOPE_HPP
