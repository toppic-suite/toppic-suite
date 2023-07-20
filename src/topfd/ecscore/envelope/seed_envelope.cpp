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

#include "common/base/mass_constant.hpp"
#include "ms/spec/peak_util.hpp"
#include "ms/env/env_base.hpp"
#include "topfd/ecscore/envelope/seed_envelope.hpp"

namespace toppic {

SeedEnvelope::SeedEnvelope(DeconvPeakPtr peak_ptr) {
  //init for parent class
  mono_mz_ = peak_ptr->getMonoMz();
  charge_ = peak_ptr->getCharge();
  double mass = peak_ptr->getMonoMass();
  EnvPtr theo_env_ptr = EnvBase::getEnvByMonoMass(mass, charge_);
  refer_idx_ = theo_env_ptr->getReferIdx();
  for (int i = 0; i < theo_env_ptr->getPeakNum(); i++) {
    EnvPeakPtr p_ptr = std::make_shared<EnvPeak>(theo_env_ptr->getMz(i),
                                                 theo_env_ptr->getInte(i));
    peak_ptr_list_.push_back(p_ptr);
  }
  // init for seed envelope
  spec_id_ = peak_ptr->getSpId();
  env_id_ = peak_ptr->getPeakId();
  seed_inte_ = peak_ptr->getIntensity();
}

SeedEnvelope::SeedEnvelope(int spec_id, int env_id, double pos, double mass,
                           double inte, int charge, std::vector<double> pos_list, 
                           std::vector<double> inte_list) {
  spec_id_ = spec_id;
  env_id_ = env_id;
  mono_mz_ = pos;
    seed_inte_ = inte;
  charge_ = charge;
  int num_peaks = pos_list.size();
  for (int i = 0; i < num_peaks; i++) {
    EnvPeakPtr p_ptr = std::make_shared<EnvPeak>(pos_list[i], inte_list[i]);
    peak_ptr_list_.push_back(p_ptr);
  }
}

SeedEnvelope::SeedEnvelope(SeedEnvelopePtr env_ptr) {
  spec_id_ = env_ptr->spec_id_;
  env_id_ = env_ptr->env_id_;
  mono_mz_ = env_ptr->mono_mz_;
  seed_inte_ = env_ptr->seed_inte_;
  charge_ = env_ptr->charge_;
  double mass = env_ptr->getMonoNeutralMass();
  EnvPtr theo_env_ptr = EnvBase::getEnvByMonoMass(mass, charge_);
  refer_idx_ = theo_env_ptr->getReferIdx();
  for (int j = 0; j < theo_env_ptr->getPeakNum(); j++) {
    EnvPeakPtr p_ptr = std::make_shared<EnvPeak>(theo_env_ptr->getMz(j),
                                                 theo_env_ptr->getInte(j));
    peak_ptr_list_.push_back(p_ptr);
  }
}

SeedEnvelope::SeedEnvelope(SeedEnvelopePtr env_ptr, int new_charge) {
  double mass = env_ptr->getMonoNeutralMass(); 
  spec_id_ = env_ptr->spec_id_;
  env_id_ = -1;
    seed_inte_ = env_ptr->seed_inte_;
  charge_ = new_charge;
  mono_mz_ = peak_util::compMz(mass, charge_);

  int old_charge = env_ptr->getCharge();
  for (auto &i: env_ptr->peak_ptr_list_) {
    double old_pos = i->getPosition();
    double neutral_mass = peak_util::compPeakNeutralMass(old_pos, old_charge);
    double new_pos = peak_util::compMz(neutral_mass, new_charge); 
    EnvPeakPtr p_ptr = std::make_shared<EnvPeak>(new_pos,
                                                 i->getIntensity());
    peak_ptr_list_.push_back(p_ptr);
  }
}

double SeedEnvelope::getReferMz() {
  double mono_mass = getMonoNeutralMass();
  double ref_mass = EnvBase::convertMonoMassToRefMass(mono_mass);
  double ref_mz = peak_util::compMz(ref_mass, charge_);
  return ref_mz;
}

void SeedEnvelope::seedRmPeaks(double min_pos, double max_pos) {
  EnvPeakPtrVec new_ptr_list;
  for (auto p: peak_ptr_list_) {
    if (p->getPosition() >= min_pos && p->getPosition() <= max_pos)
      new_ptr_list.push_back(p);
  }
    peak_ptr_list_ = new_ptr_list;
}


void SeedEnvelope::seedShiftIsotope(double shift_num) {
  double shift_mass = shift_num * mass_constant::getIsotopeMass();
  double shift_mz = shift_mass / charge_;
  mono_mz_ = mono_mz_ + shift_mz;
  for (auto &p: peak_ptr_list_)
    p->setPosition(p->getPosition() + shift_mz);
}

std::string SeedEnvelope::getString() {
  std::string header = "Spec ID: " + std::to_string(spec_id_) + " " +
                       "Env ID: " + std::to_string(env_id_) + " " +
                       "Pos: " + std::to_string(mono_mz_) + " " +
                       "Inte: " + std::to_string(seed_inte_) + " " +
                       "Charge: " + std::to_string(charge_) + "\n";
  std::string peaks = "(";
  for (auto peak: peak_ptr_list_)
    peaks = peaks + "(" + std::to_string(peak->getPosition()) + ", " 
      + std::to_string(peak->getIntensity()) + "), ";
  peaks = peaks + ")\n";
  return header + peaks;
}

}
