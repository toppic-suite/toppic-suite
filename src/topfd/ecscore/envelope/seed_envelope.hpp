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

#ifndef TOPPIC_TOPFD_ECSCORE_ENVELOPE_SEED_ENVELOPE_HPP
#define TOPPIC_TOPFD_ECSCORE_ENVELOPE_SEED_ENVELOPE_HPP

#include <vector>
#include <algorithm>

#include "ms/spec/deconv_ms.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/env_base.hpp"
//#include "topfd/ecscore/envelope/env_simple_peak.hpp"

namespace toppic {

class SeedEnvelope;

typedef std::shared_ptr<SeedEnvelope> SeedEnvelopePtr;

class SeedEnvelope : public Env {
 public:

  int getSpecId() const { return spec_id_; }

  void setSpecId(int spec_id) { spec_id_ = spec_id; }

  int getEnvId() const { return env_id_; }

  double getSeedInte() const {return seed_inte_;}

  static bool cmpSeedInteDec(const SeedEnvelopePtr a, const SeedEnvelopePtr b) {
    return a->getSeedInte() > b->getSeedInte(); }

  std::string getString();

  // to review
  SeedEnvelope(DeconvPeakPtr peak_ptr);

  SeedEnvelope(int spec_id, int env_id, double pos, double mass, double inte, int charge,
               std::vector<double> pos_list, std::vector<double> inte_list);

  SeedEnvelope(SeedEnvelopePtr env_ptr);

  SeedEnvelope(SeedEnvelopePtr env_ptr, int new_charge);

  void seedRmPeaks(double min_pos, double max_pos);

  void seedShiftIsotope(double shift_num);

  double getReferMz();

 private:
  int spec_id_;
  int env_id_;
  double seed_inte_;
};

typedef std::vector<SeedEnvelopePtr> SeedEnvelopePtrVec;
typedef std::vector<SeedEnvelopePtrVec> SeedEnvelopePtr2D;

}

#endif //TOPPIC_SEED_ENVELOPE_HPP
