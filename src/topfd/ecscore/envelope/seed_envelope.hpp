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
#include "topfd/ecscore/spectrum/env_simple_peak.hpp"

namespace toppic {

class SeedEnvelope;

typedef std::shared_ptr<SeedEnvelope> SeedEnvelopePtr;

class SeedEnvelope {
 public:
  SeedEnvelope();

  SeedEnvelope(DeconvPeakPtr peak_ptr);

  SeedEnvelope(int spec_id, int env_id, double pos, double mass, double inte, int charge,
               std::vector<double> pos_list, std::vector<double> inte_list);

  SeedEnvelope(SeedEnvelopePtr env_ptr);

  SeedEnvelope(MsHeaderPtr header_ptr);

  std::vector<double> getPosList();

  std::vector<double> getInteList();

  void rmPeaks(double min_pos, double max_pos);

  void changeCharge(int new_charge);

  void shift(double shift_num);

  double getMaxPos() { return peak_ptr_list_[peak_ptr_list_.size() - 1]->getPosition(); }

  int getPeakNum() { return peak_ptr_list_.size(); }

  EnvSimplePeakPtr get_peak_ptr(int idx) { return peak_ptr_list_[idx]; }

  int getSpecId() const { return spec_id_; }

  void setSpecId(int spec_id) { spec_id_ = spec_id; }

  int getEnvId() const { return env_id_; }

  void setEnvId(int env_id) { env_id_ = env_id; }

  double getPos() const { return pos_; }

  void setPos(double pos) { pos_ = pos; }

  double getMass() const { return mass_; }

  void setMass(double mass) { mass_ = mass; }

  double getInte() const { return inte_; }

  void setInte(double inte) { inte_ = inte; }

  int getCharge() const { return charge_; }

  void setCharge(int charge) { charge_ = charge; }

  const EnvSimplePeakPtrVec &getPeakList() const { return peak_ptr_list_; }

  void setPeakList(EnvSimplePeakPtrVec peak_ptr_list); 

  std::string getString();

  static bool cmpInteDec(const SeedEnvelopePtr a, const SeedEnvelopePtr b) { 
    return a->getInte() > b->getInte(); }

 private:
  int spec_id_;
  int env_id_;
  double pos_;
  double mass_;
  double inte_;
  int charge_;
  EnvSimplePeakPtrVec peak_ptr_list_;
};

typedef std::vector<SeedEnvelopePtr> SeedEnvelopePtrVec;

}

#endif //TOPPIC_SEED_ENVELOPE_HPP
