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

#include "env_collection.hpp"

namespace toppic {

  EnvCollection::EnvCollection() {
    seed_env_ = SeedEnvelope();
    min_charge_ = -1;
    max_charge_ = -1;
    start_spec_id_ = -1;
    end_spec_id_ = -1;
  }

  EnvCollection::EnvCollection(const EnvCollection &ec) {
    seed_env_ = ec.seed_env_;
    min_charge_ = ec.min_charge_;
    max_charge_ = ec.max_charge_;
    start_spec_id_ = ec.start_spec_id_;
    end_spec_id_ = ec.end_spec_id_;
    ecscore_ = ec.ecscore_;
    for (auto &i: ec.env_set_list_)
      env_set_list_.push_back(i);
    for (auto &i: ec.exp_inte_sum_list_)
      exp_inte_sum_list_.push_back(i);
  }

  EnvCollection::EnvCollection(const SeedEnvelope &env, const std::vector<EnvSet> &env_set_list,
                               int min_charge, int max_charge, int start_spec_id, int end_spec_id) {
    seed_env_ = env;
    for (auto &i: env_set_list) env_set_list_.push_back(i);
    min_charge_ = min_charge;
    max_charge_ = max_charge;
    start_spec_id_ = start_spec_id;
    end_spec_id_ = end_spec_id;
    exp_inte_sum_list_ = comp_exp_inte_sum_list();
  }

  std::vector<double> EnvCollection::comp_exp_inte_sum_list() {
    int peak_num = seed_env_.get_peak_num();
    std::vector<double> sum_list(peak_num, 0);
    for (auto &env_set: env_set_list_) {
      if (env_set.isEmpty())
        continue;
      std::vector<double> cur_sum_list = env_set.comp_exp_inte_sum_list();
      if (peak_num != static_cast<int>(cur_sum_list.size()))
        std::cout << "peak number " << peak_num << " cur list len " << cur_sum_list.size();
      for (int i = 0; i < peak_num; i++)
        sum_list[i] = sum_list[i] + cur_sum_list[i];
    }
    return sum_list;
  }

  bool EnvCollection::isEmpty() {
    if (seed_env_.isEmpty() && env_set_list_.empty() && min_charge_ == -1 && max_charge_ == -1 &&
        start_spec_id_ == -1 && end_spec_id_ == -1)
      return true;
    return false;
  }

  std::vector<int> EnvCollection::getChargeList() {
    std::vector<int> charge_list;
    for (auto es: env_set_list_)
      charge_list.push_back(es.getCharge());
    return charge_list;
  }

  void EnvCollection::refine_mono_mass() {
    double weight = 0;
    double weight_mz_error = 0;
    for (auto &env_set: env_set_list_) {
      double cur_weight = 0, cur_weight_mz_error = 0;
      env_set.get_weight_mz_error(&cur_weight, &cur_weight_mz_error);
      weight = weight + cur_weight;
      weight_mz_error = weight_mz_error + cur_weight_mz_error;
    }
    if (weight > 0) {
      double mz_error = weight_mz_error / weight;
      seed_env_.shift(mz_error * seed_env_.getCharge());
    }
    else {
      std::cout << "ERROR 0 weight in refine_mono_mass" << std::endl;
    }
  }

  double EnvCollection::get_intensity(double snr, double noise_inte) {
    double inte = 0;
    for (auto env_set: env_set_list_) {
      double tmp_inte = env_set.comp_intensity(snr, noise_inte);
      inte = inte + tmp_inte;
    }
    return inte;
  }

  double EnvCollection::get_apex_elution_time(spec_list &spectra_list) {
    int ref_spec_id = seed_env_.getSpecId();
    return spectra_list[ref_spec_id].getRt();
  }

  double EnvCollection::get_elution_length(spec_list &spectra_list) {
    return get_max_elution_time(spectra_list) - get_min_elution_time(spectra_list);
  }

  void EnvCollection::remove_peak_data(PeakMatrix &peak_matrix) {
    for (auto env_set: env_set_list_) {
      if (!env_set.isEmpty())
        env_set.remove_peak_data(peak_matrix);
    }
  }

  std::vector<std::vector<double>> EnvCollection::get_seed_theo_map(PeakMatrix &peak_matrix, double snr) {
    double noise_inte = peak_matrix.get_min_inte();
    EnvSet env_set = get_seed_env_set();
    std::vector<std::vector<double>> map = env_set.get_map(snr, noise_inte);
    return map;
  }

  toppic::EnvSet EnvCollection::get_seed_env_set() {
    EnvSet env_set = EnvSet();
    for (const auto &es: env_set_list_) {
      SeedEnvelope es_seed_env = es.getSeedEnv();
      if (es_seed_env.getCharge() == seed_env_.getCharge())
        env_set = EnvSet(es);
    }
    return env_set;
  }

  void EnvCollection::setEnvSetList(const std::vector<EnvSet> &env_set_list) {
    for (auto &i: env_set_list)
      env_set_list_.push_back(i);
  }

  void EnvCollection::setExpInteSumList(const std::vector<double> &expInteSumList) {
    exp_inte_sum_list_.clear();
    for (auto &p: expInteSumList)
      exp_inte_sum_list_.push_back(p);
  }
}