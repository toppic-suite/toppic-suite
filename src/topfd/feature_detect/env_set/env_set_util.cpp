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

#include "env_set_util.hpp"

namespace toppic {
namespace env_set_util {
    ExpPeak _pick_exp_peak(PeakMatrix& peak_matrix, SimplePeak& seed_peak, int sp_id, double mass_tol) {
      // get peaks within mass tolerance
      ExpPeak result_peak = ExpPeak();
      double max_inte = std::numeric_limits<double>::min();
      double pos = seed_peak.getPos();
      for (int idx = seed_peak.getStartIdx(); idx < seed_peak.getEndIdx() + 1; idx++) {
        std::vector<ExpPeak> bin_peaks = peak_matrix.get_bin_peak(sp_id, idx);
        for (const auto& matrix_peak : bin_peaks) {
          double mass_diff = std::abs(pos - matrix_peak.getPos());
          if ( mass_diff < mass_tol && matrix_peak.getInte() > max_inte) {
            result_peak = matrix_peak; ////////////////////// this is slow!!!!
            max_inte = matrix_peak.getInte();
          }
        }
      }
      return result_peak;
    }

    ExpEnvelope get_match_exp_env(PeakMatrix &peak_matrix, SeedEnvelope &seed_env, int sp_id, double mass_tol) {
      std::vector<ExpPeak> peak_list;
      std::vector<SimplePeak> peaks = seed_env.getPeakList();
      ExpPeak peak = ExpPeak();
      for (auto& seed_peak : peaks) {
        if (seed_peak.getStartIdx() > -1 && seed_peak.getEndIdx() > -1)
          peak = _pick_exp_peak(peak_matrix, seed_peak, sp_id, mass_tol);
        peak_list.push_back(peak);
      }
      ExpEnvelope exp_env = ExpEnvelope(sp_id, peak_list);
      return exp_env;
    }

  bool check_valid_env_set_seed_env(PeakMatrix& peak_matrix, EnvSet& env_set, int max_miss_peak) {
    std::vector<double> theo_envelope_inte = env_set.get_theo_distribution_inte();
    int refer_idx = std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) - theo_envelope_inte.begin();
    int base_idx = env_set.getBaseSpecId();
    int start_idx = std::max(base_idx-1, 0);
    int end_idx = std::min(base_idx+1, peak_matrix.get_spec_num() -1);
    std::vector<ExpEnvelope> env_list = env_set.getExpEnvList();
    bool valid = true;
    for (auto &exp_env : env_list) {
      if (exp_env.getSpecId() >= start_idx and exp_env.getSpecId() <= end_idx)
        if (exp_env.get_match_peak_num(refer_idx) < max_miss_peak)
          valid = false;
    }
    return valid;
  }

  bool check_valid_env_set(PeakMatrix& peak_matrix, EnvSet& env_set) {
    bool valid = true;
    int elems = 0;
    std::vector<double> env_xic = env_set.getXicEnvIntes();
    for (double inte : env_xic)
      if (inte > 0) elems++;
    if (elems < 2) valid = false;
    return valid;
  }

  void remove_non_match_envs(std::vector<ExpEnvelope> &env_list, int refer_idx) {
    int idx = env_list.size() - 1;
    while (idx >= 0) {
      ExpEnvelope env = env_list[idx];
      if (env.get_match_peak_num(refer_idx) < 2)
        env_list.erase(env_list.begin() + idx);
      else
        return;
      idx = idx - 1;
    }
  }

  void comp_peak_start_end_idx(PeakMatrix &peak_matrix, SeedEnvelope &seed_env, double error_tole) {
    std::vector<SimplePeak> peak_list = seed_env.getPeakList();
    for (auto &peak : peak_list) {
      double mz = peak.getPos();
      int start_idx = peak_matrix.get_index(mz - error_tole);
      if (start_idx < 0)
        start_idx = 0;
      int end_idx = peak_matrix.get_index(mz + error_tole);
      if (end_idx >= peak_matrix.get_bin_num())
        end_idx = peak_matrix.get_bin_num() - 1;
      peak.setStartIdx(start_idx);
      peak.setEndIdx(end_idx);
    }
    seed_env.setPeakList(peak_list);
  }

  EnvSet find_env_set(PeakMatrix &peak_matrix, SeedEnvelope &env, int start_spec_id, int end_spec_id, FeatureParaPtr para_ptr, double sn_ratio) {
    double noise_inte_level = peak_matrix.get_min_inte();

    ExpEnvelope empty_exp_env = ExpEnvelope();
    std::vector<double> theo_envelope_inte = env.get_inte_list();
    int num_theo_env_peaks = theo_envelope_inte.size();
    int refer_idx = std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) - theo_envelope_inte.begin();
    int base_idx = env.getSpecId();
    int miss_num = 0;

    std::vector<ExpEnvelope> back_env_list;
    for (int idx = base_idx; idx >= start_spec_id; idx--) {
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix, env, idx, para_ptr->mass_tole_);
      std::vector<double> experimental_envelope_inte = exp_env.get_inte_list();
      double inte_ratio = env_utils::calcInteRatio_scan(theo_envelope_inte, experimental_envelope_inte);
      for (int i = 0; i < num_theo_env_peaks; i++) {
        double peak_inte = theo_envelope_inte[i];
        if ((inte_ratio * peak_inte) < (noise_inte_level * sn_ratio))
          exp_env.setExpEnvListPeak(ExpPeak(), i);
      }
      back_env_list.push_back(exp_env);
      if (exp_env.get_match_peak_num(refer_idx) < para_ptr->max_miss_peak_)
        miss_num = miss_num + 1;
      else
        miss_num = 0;
      if (miss_num >=  para_ptr->max_miss_env_)
        break;
    }
    env_set_util::remove_non_match_envs(back_env_list, refer_idx);

    std::vector<ExpEnvelope> forw_env_list;
    for (int idx = base_idx + 1; idx <= end_spec_id; idx++) {
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix, env, idx,  para_ptr->mass_tole_);
      std::vector<double> experimental_envelope_inte = exp_env.get_inte_list();
      double inte_ratio = env_utils::calcInteRatio_scan(theo_envelope_inte, experimental_envelope_inte);
      for (int i = 0; i < num_theo_env_peaks; i++) {
        double peak_inte = theo_envelope_inte[i];
        if ((inte_ratio * peak_inte) < (noise_inte_level * sn_ratio))
          exp_env.setExpEnvListPeak(ExpPeak(), i);
      }
      forw_env_list.push_back(exp_env);
      if (exp_env.get_match_peak_num(refer_idx) <  para_ptr->max_miss_peak_)
        miss_num = miss_num + 1;
      else
        miss_num = 0;
      if (miss_num >= para_ptr->max_miss_env_)
        break;
    }
    env_set_util::remove_non_match_envs(forw_env_list, refer_idx);
    // merge
    std::reverse(back_env_list.begin(), back_env_list.end());
    back_env_list.insert(back_env_list.end(), forw_env_list.begin(), forw_env_list.end());
    if (back_env_list.empty()) return EnvSet();
    start_spec_id = back_env_list[0].getSpecId();
    end_spec_id = back_env_list[back_env_list.size() - 1].getSpecId();
    if ((end_spec_id - start_spec_id + 1) < 2) return EnvSet();
    EnvSet env_set = EnvSet(env, back_env_list, start_spec_id, end_spec_id, noise_inte_level, sn_ratio);
    return env_set;
  }

  EnvSet get_env_set(PeakMatrix& peak_matrix, SeedEnvelope& env, FeatureParaPtr para_ptr, double sn_ratio) {
    double noise_inte_level = peak_matrix.get_min_inte();
    std::vector<double> theo_envelope_inte = env.get_inte_list();
    int num_theo_env_peaks = theo_envelope_inte.size();
    int refer_idx = std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) - theo_envelope_inte.begin();

    // search backward
    std::vector<ExpEnvelope> back_env_list;
    int idx = env.getSpecId();
    int miss_num = 0;
    while (idx >= 0) {
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix, env, idx, para_ptr->mass_tole_);
      std::vector<double> experimental_envelope_inte = exp_env.get_inte_list();
      double inte_ratio = env_utils::calcInteRatio_scan(theo_envelope_inte, experimental_envelope_inte);
      for (int i = 0; i < num_theo_env_peaks; i++) {
        double peak_inte = theo_envelope_inte[i];
        if ((inte_ratio * peak_inte) < (noise_inte_level * sn_ratio))
          exp_env.setExpEnvListPeak(ExpPeak(), i);
      }
      back_env_list.push_back(exp_env);
      if (exp_env.get_match_peak_num(refer_idx) < para_ptr->max_miss_peak_)
        miss_num = miss_num + 1;
      else
        miss_num = 0;
      if (miss_num >= para_ptr->max_miss_env_)
        break;
      idx = idx - 1;
    }
    env_set_util::remove_non_match_envs(back_env_list, refer_idx);

    // search forward
    std::vector<ExpEnvelope> forw_env_list;
    idx = env.getSpecId() + 1;
    miss_num = 0;
    while (idx < peak_matrix.get_spec_num()) {
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix, env, idx, para_ptr->mass_tole_);
      std::vector<double> experimental_envelope_inte = exp_env.get_inte_list();
      double inte_ratio = env_utils::calcInteRatio_scan(theo_envelope_inte, experimental_envelope_inte);
      for (int i = 0; i < num_theo_env_peaks; i++) {
        double peak_inte = theo_envelope_inte[i];
        if ((inte_ratio * peak_inte) < (noise_inte_level * sn_ratio))
          exp_env.setExpEnvListPeak(ExpPeak(), i);
      }
      forw_env_list.push_back(exp_env);
      if (exp_env.get_match_peak_num(refer_idx) < para_ptr->max_miss_peak_)
        miss_num = miss_num + 1;
      else
        miss_num = 0;
      if (miss_num >= para_ptr->max_miss_env_)
        break;
      idx = idx + 1;
    }
    env_set_util::remove_non_match_envs(forw_env_list, refer_idx);
    // merge results
    std::reverse(back_env_list.begin(), back_env_list.end());
    back_env_list.insert(back_env_list.end(), forw_env_list.begin(), forw_env_list.end());
    if (back_env_list.empty()) return EnvSet();
    int start_spec_id = back_env_list[0].getSpecId();
    int end_spec_id = back_env_list[back_env_list.size() - 1].getSpecId();
    if ((end_spec_id - start_spec_id) < para_ptr->match_peak_tole_) return EnvSet();
    EnvSet env_set = EnvSet(env, back_env_list, start_spec_id, end_spec_id, noise_inte_level, sn_ratio);
    return env_set;
  }
}
}