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
#include "topfd/ecscore/spectrum/env_simple_peak.hpp"
#include "topfd/ecscore/envelope/seed_envelope.hpp"

namespace toppic {

EnvelopePtr getTheoEnv(double mono_mass, int charge) {
  EnvelopePtr ref_env_ptr = toppic::EnvBase::getStaticEnvByMonoMass(mono_mass);
  if (ref_env_ptr == nullptr) {
    return nullptr;
  }
  double mono_mz = peak_util::compMz(mono_mass, charge);
  EnvelopePtr theo_env_ptr = ref_env_ptr->distrToTheoMono(mono_mz, charge);
  return theo_env_ptr;
}

SeedEnvelope::SeedEnvelope(DeconvPeakPtr peak_ptr) {
  spec_id_ = peak_ptr->getSpId();
  env_id_ = peak_ptr->getId();
  mass_ = peak_ptr->getMonoMass();
  pos_ = peak_ptr->getMonoMz();
  inte_ = peak_ptr->getIntensity();
  charge_ = peak_ptr->getCharge();
  EnvelopePtr theo_env_ptr = getTheoEnv(mass_, charge_);
  for (int j = 0; j < theo_env_ptr->getPeakNum(); j++) {
    EnvSimplePeakPtr p_ptr = std::make_shared<EnvSimplePeak>(theo_env_ptr->getMz(j), 
                                                           theo_env_ptr->getIntensity(j));
    peak_ptr_list_.push_back(p_ptr);
  }
}

SeedEnvelope::SeedEnvelope(MsHeaderPtr header_ptr) {
  spec_id_ = header_ptr->getMsOneId();
  env_id_ = -1;
  mass_ = header_ptr->getPrecMonoMass();
  pos_ = header_ptr->getPrecMonoMz();
  inte_ = header_ptr->getPrecInte();
  charge_ = header_ptr->getPrecCharge();
  EnvelopePtr theo_env_ptr = getTheoEnv(mass_, charge_);
  for (int j = 0; j < theo_env_ptr->getPeakNum(); j++) {
    EnvSimplePeakPtr p_ptr = std::make_shared<EnvSimplePeak>(theo_env_ptr->getMz(j), 
                                                             theo_env_ptr->getIntensity(j));
    peak_ptr_list_.push_back(p_ptr);
  }
}

SeedEnvelope::SeedEnvelope(int spec_id, int env_id, double pos, double mass, double inte, int charge,
                           std::vector<double> pos_list, std::vector<double> inte_list) {
  spec_id_ = spec_id;
  env_id_ = env_id;
  pos_ = pos;
  mass_ = mass;
  inte_ = inte;
  charge_ = charge;
  int num_peaks = pos_list.size();
  for (int i = 0; i < num_peaks; i++) {
    EnvSimplePeakPtr p_ptr = std::make_shared<EnvSimplePeak>(pos_list[i], inte_list[i]); 
    peak_ptr_list_.push_back(p_ptr);
  }
}

SeedEnvelope::SeedEnvelope(SeedEnvelopePtr env_ptr) {
  spec_id_ = env_ptr->spec_id_;
  env_id_ = env_ptr->env_id_;
  pos_ = env_ptr->pos_;
  mass_ = env_ptr->mass_;
  inte_ = env_ptr->inte_;
  charge_ = env_ptr->charge_;
  EnvelopePtr theo_env_ptr = getTheoEnv(mass_, charge_);
  for (int j = 0; j < theo_env_ptr->getPeakNum(); j++) {
    EnvSimplePeakPtr p_ptr = std::make_shared<EnvSimplePeak>(theo_env_ptr->getMz(j), 
                                                           theo_env_ptr->getIntensity(j));
    peak_ptr_list_.push_back(p_ptr);
  }
}

SeedEnvelope::SeedEnvelope(SeedEnvelopePtr env_ptr, int new_charge) {
  spec_id_ = env_ptr->spec_id_;
  env_id_ = -1; 
  mass_ = env_ptr->mass_;
  inte_ = env_ptr->inte_;
  charge_ = new_charge;
  pos_ = peak_util::compMz(mass_, charge_);

  int old_charge = env_ptr->getCharge();
  for (auto &i: env_ptr->peak_ptr_list_) {
    double old_pos = i->getPosition();
    double neutral_mass = peak_util::compPeakNeutralMass(old_pos, old_charge);
    double new_pos = peak_util::compMz(neutral_mass, new_charge); 
    EnvSimplePeakPtr p_ptr = std::make_shared<EnvSimplePeak>(new_pos,  
                                                             i->getIntensity()); 
    peak_ptr_list_.push_back(p_ptr);
  }
}


std::vector<double> SeedEnvelope::getPosList() {
  std::vector<double> pos_list;
  for (auto p: peak_ptr_list_)
    pos_list.push_back(p->getPosition());
  return pos_list;
}

std::vector<double> SeedEnvelope::getInteList() {
  std::vector<double> inte_list;
  for (auto p: peak_ptr_list_)
    inte_list.push_back(p->getIntensity());
  return inte_list;
}

void SeedEnvelope::rmPeaks(double min_pos, double max_pos) {
  EnvSimplePeakPtrVec new_ptr_list;
  for (auto p: peak_ptr_list_) {
    if (p->getPosition() >= min_pos && p->getPosition() <= max_pos)
      new_ptr_list.push_back(p);
  }
  peak_ptr_list_ = new_ptr_list;
}

void SeedEnvelope::changeCharge(int new_charge) {
  for (auto &p: peak_ptr_list_) {
    double new_pos = (((p->getPosition() - mass_constant::getProtonMass()) * charge_) / new_charge) 
      + mass_constant::getProtonMass();
    p->setPosition(new_pos);
  }
  pos_ = (((pos_ - mass_constant::getProtonMass()) * charge_) / new_charge) 
    + mass_constant::getProtonMass();
  charge_ = new_charge;
}

void SeedEnvelope::shift(double shift_num) {
  double shift_mass = shift_num * mass_constant::getIsotopeMass();
  double shift_mz = shift_mass / charge_;
  pos_ = pos_ + shift_mz;
  mass_ = mass_ + shift_mass;
  for (auto &p: peak_ptr_list_)
    p->setPosition(p->getPosition() + shift_mz);
}

void SeedEnvelope::setPeakList(EnvSimplePeakPtrVec new_peak_list) {
  peak_ptr_list_ = new_peak_list;
}

std::string SeedEnvelope::getString() {
  std::string header = "Spec ID: " + std::to_string(spec_id_) + " " +
    "Env ID: " + std::to_string(env_id_) + " " +
    "Pos: " + std::to_string(pos_) + " " +
    "Inte: " + std::to_string(inte_) + " " +
    "Mass: " + std::to_string(mass_) + " " +
    "Charge: " + std::to_string(charge_) + "\n";
  std::string peaks = "(";
  for (auto peak: peak_ptr_list_)
    peaks = peaks + "(" + std::to_string(peak->getPosition()) + ", " 
      + std::to_string(peak->getIntensity()) + "), ";
  peaks = peaks + ")\n";
  return header + peaks;
}

}
