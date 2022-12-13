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

#include "seed_envelope.hpp"

namespace toppic {
  std::shared_ptr<toppic::Envelope> _get_env(double mono_mass, int charge, double mono_mz) {
    toppic::EnvelopePtr ref_env_ptr = toppic::EnvBase::getStaticEnvByMonoMass(mono_mass);
    if (ref_env_ptr == nullptr) {
      return nullptr;
    }
    toppic::EnvelopePtr theo_env_ptr = ref_env_ptr->distrToTheoMono(mono_mz, charge);
    return theo_env_ptr;
  }

  SeedEnvelope::SeedEnvelope() {
    spec_id_ = -1;
    env_id_ = -1;
    pos_ = -1;
    mass_ = -1;
    inte_ = -1;
    charge_ = -1;
  }

  SeedEnvelope::SeedEnvelope(DeconvPeakPtr &p) {
    spec_id_ = p->getSpId();
    env_id_ = p->getId();
    mass_ = p->getMonoMass();
    pos_ = p->getMonoMz();
    inte_ = p->getIntensity();
    charge_ = p->getCharge();
    std::shared_ptr<toppic::Envelope> theo_env = _get_env(mass_, charge_, pos_);
    for (int j = 0; j < theo_env->getPeakNum(); j++)
      peak_list_.push_back(SimplePeak(theo_env->getMz(j), theo_env->getIntensity(j)));
  }

  SeedEnvelope::SeedEnvelope(MsHeaderPtr &hh) {
    spec_id_ = hh->getMsOneId();
    env_id_ = 0;
    mass_ = hh->getPrecMonoMass();
    pos_ = hh->getPrecMonoMz();
    inte_ = hh->getPrecInte();
    charge_ = hh->getPrecCharge();
    std::shared_ptr<toppic::Envelope> theo_env = _get_env(mass_, charge_, pos_);
    for (int j = 0; j < theo_env->getPeakNum(); j++)
      peak_list_.push_back(SimplePeak(theo_env->getMz(j), theo_env->getIntensity(j)));
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
    for (int i = 0; i < num_peaks; i++)
      peak_list_.push_back(SimplePeak(pos_list[i], inte_list[i]));
  }

  SeedEnvelope::SeedEnvelope(const toppic::SeedEnvelope &env) {
    spec_id_ = env.spec_id_;
    env_id_ = env.env_id_;
    pos_ = env.pos_;
    mass_ = env.mass_;
    inte_ = env.inte_;
    charge_ = env.charge_;
    for (const auto &i: env.peak_list_)
      peak_list_.push_back(SimplePeak(i.getPos(), i.getInte()));
  }

  std::vector<double> SeedEnvelope::get_pos_list() {
    std::vector<double> pos_list;
    for (auto p: peak_list_)
      pos_list.push_back(p.getPos());
    return pos_list;
  }

  std::vector<double> SeedEnvelope::get_inte_list() {
    std::vector<double> inte_list;
    for (auto p: peak_list_)
      inte_list.push_back(p.getInte());
    return inte_list;
  }

  void SeedEnvelope::rm_peaks(double min_pos, double max_pos) {
    std::vector<SimplePeak> peak_list;
    for (auto p: peak_list_) {
      if (p.getPos() >= min_pos and p.getPos() <= max_pos)
        peak_list.push_back(p);
      peak_list_ = peak_list;
    }
  }

  void SeedEnvelope::change_charge(int new_charge) {
    for (auto &p: peak_list_) {
      double new_pos = (((p.getPos() - get_proton_mass()) * charge_) / new_charge) + get_proton_mass();
      p.setPos(new_pos);
    }
    pos_ = (((pos_ - get_proton_mass()) * charge_) / new_charge) + get_proton_mass();
    charge_ = new_charge;
  }

  SeedEnvelope SeedEnvelope::get_new_charge_env(int new_charge) {
    SeedEnvelope new_env = *this;
    new_env.change_charge(new_charge);
    return new_env;
  }

  void SeedEnvelope::shift(double shift_num) {
    double shift_mass = shift_num * get_isotope_mass();
    double shift_mz = shift_mass / charge_;
    pos_ = pos_ + shift_mz;
    mass_ = mass_ + shift_mass;
    for (auto &p: peak_list_)
      p.setPos(p.getPos() + shift_mz);
  }

  void SeedEnvelope::setPeakList(const std::vector<SimplePeak> &peakList) {
    peak_list_.clear();
    for (auto &peak: peakList)
      peak_list_.push_back(peak);
  }

  bool SeedEnvelope::isEmpty() {
    if (spec_id_ == -1 && env_id_ == -1 && pos_ == -1 && mass_ == -1 && inte_ == -1 && charge_ == -1 &&
        peak_list_.empty())
      return true;
    return false;
  }

  std::string SeedEnvelope::getString() {
    std::string header = "Spec ID: " + std::to_string(spec_id_) + " " +
                         "Env ID: " + std::to_string(env_id_) + " " +
                         "Pos: " + std::to_string(pos_) + " " +
                         "Inte: " + std::to_string(inte_) + " " +
                         "Mass: " + std::to_string(mass_) + " " +
                         "Charge: " + std::to_string(charge_) + "\n";
    std::string peaks = "(";
    for (auto peak: peak_list_)
      peaks = peaks + "(" + std::to_string(peak.getPos()) + ", " + std::to_string(peak.getInte()) + "), ";
    peaks = peaks + ")\n";
    return header + peaks;
  }
}