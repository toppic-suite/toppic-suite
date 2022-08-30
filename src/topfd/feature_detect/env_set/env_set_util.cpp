#include <valarray>
#include <iostream>
#include "env_set_util.hpp"

namespace toppic {
namespace env_set_util {
  ExpPeak pick_exp_peak(PeakRow row, SimplePeak seed_peak, double mass_tol) {
    // get peaks within mass tolerance
    ExpPeak result_peak = ExpPeak();
    double max_inte = -1000000;
    double pos = seed_peak.getPos();
    for (int idx = seed_peak.getStartIdx(); idx < seed_peak.getEndIdx() + 1; idx++) {
      std::vector<std::vector<ExpPeak>> rows = row.getRow();
      for (auto matrix_peak : rows[idx]) {
        if (std::abs(pos - matrix_peak.getPos()) < mass_tol && (matrix_peak.getInte() > max_inte)) {
          result_peak = ExpPeak(matrix_peak);
          max_inte = matrix_peak.getInte();
        }
      }
    }
    return result_peak;
  }

  ExpEnvelope get_match_exp_env(PeakRow peak_row, SeedEnvelope seed_env, double mass_tol) {
    std::vector<ExpPeak> peak_list;
    std::vector<SimplePeak> peaks = seed_env.getPeakList();
    for (auto seed_peak : peaks) {
      ExpPeak peak = pick_exp_peak(peak_row, seed_peak, mass_tol);
      peak_list.push_back(peak);
    }
    ExpEnvelope exp_env = ExpEnvelope(peak_row.getSpecID(), peak_list);
    return exp_env;
  }

  SeedEnvelope preprocess_env(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tol, bool *valid) {
    SeedEnvelope env(seed_env);
    double min_mz = peak_matrix.get_min_mz() - mass_tol;
    double max_mz = peak_matrix.get_max_mz() + mass_tol;
    env.rm_peaks(min_mz, max_mz);
    // if peak number is low, stop
    *valid = true;
    if (env.get_peak_num() <= 2)
      *valid = false;
    if (env.getSpecId() >= peak_matrix.get_spec_num()) {
      std::cout << "spec id " << seed_env.getSpecId() << " is out of range!" << std::endl;
      *valid = false;
    }
    return env;
  }

    bool check_valid_env_set(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tol, double match_peak_num_tole) {
    bool valid = false;
    SeedEnvelope env = env_set_util::preprocess_env(peak_matrix, seed_env, mass_tol, &valid);
    if (!valid)
      return false;
    env.keep_top_three();
    env_set_util::comp_peak_start_end_idx(peak_matrix, env.getPeakList(), mass_tol);
    int base_idx = env.getSpecId();
    int start_idx = std::max(base_idx - 1, 0);
    int end_idx = std::min(base_idx + 1, peak_matrix.get_spec_num() - 1);

    valid = true;
    std::vector<ExpEnvelope> env_list;
    for (int idx = start_idx; idx < end_idx + 1; idx++) {
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix.get_row(idx), env, mass_tol);
      env_list.push_back(exp_env);
      if (exp_env.get_match_peak_num() < match_peak_num_tole)
        valid = false;
    }
    return valid;
  }

  void remove_non_match_envs(std::vector<ExpEnvelope> env_list) {
    int idx = env_list.size() - 1;
    while (idx >= 0) {
      ExpEnvelope env = env_list[idx];
      if (env.get_match_peak_num() < 2)
        env_list.erase(env_list.begin() + idx);
        // env_list.pop();
      else
        return;
      idx = idx - 1;
    }
  }

  void comp_peak_start_end_idx(PeakMatrix peak_matrix, std::vector<SimplePeak> peak_list, double error_tole) {
    for (auto peak : peak_list) {
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
  }



  EnvSet get_env_set_by_three_peaks(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tol, double max_miss_env) {
    bool valid;
    SeedEnvelope env = env_set_util::preprocess_env(peak_matrix, seed_env, mass_tol, &valid);
    if (!valid)
      return EnvSet();
    env.keep_top_three();
    comp_peak_start_end_idx(peak_matrix, env.getPeakList(), mass_tol);

    // search backward
    std::vector<ExpEnvelope> back_env_list;
    int idx = env.getSpecId();
    int miss_num = 0;
    while (idx >= 0) {
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix.get_row(idx), env, mass_tol);
      back_env_list.push_back(exp_env);
      if (exp_env.get_match_peak_num() < 2)
        miss_num = miss_num + 1;
      else
        miss_num = 0;
      if (miss_num >= max_miss_env)
        break;
      idx = idx - 1;
    }
    env_set_util::remove_non_match_envs(back_env_list);

    // search forward
    std::vector<ExpEnvelope> forw_env_list;
    idx = env.getSpecId() + 1;
    miss_num = 0;
    while (idx < peak_matrix.get_spec_num()) {
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix.get_row(idx), env, mass_tol);
      forw_env_list.push_back(exp_env);
      if (exp_env.get_match_peak_num() < 2)
        miss_num = miss_num + 1;
      else
        miss_num = 0;
      if (miss_num >= max_miss_env)
        break;
      idx = idx + 1;
    }
    env_set_util::remove_non_match_envs(forw_env_list);

    // merge results
    std::reverse(back_env_list.begin(), back_env_list.end());
    back_env_list.insert(back_env_list.end(), forw_env_list.begin(), forw_env_list.end());
    if (back_env_list.empty())
      return EnvSet();
    int start_spec_id = back_env_list[0].getSpecId();
    int end_spec_id = back_env_list[back_env_list.size() - 1].getSpecId();
    EnvSet env_set = EnvSet(env, back_env_list, start_spec_id, end_spec_id);
    return env_set;
  }

  void print_env(ExpEnvelope exp_env, double ratio) {
    std::cout << exp_env.getSpecId() << " " << std::endl;
    std::vector<ExpPeak> peak_list = exp_env.getExpEnvList();
    for (auto p : peak_list) {
      if (!p.isEmpty())
        std::cout << p.getPos() << " " << p.getInte() * ratio << std::endl;
      else
        std::cout << "NONE , NONE" << std::endl;
    }
    std::cout << std::endl;
  }

  EnvSet find_env_set(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tol, int start_spec_id, int end_spec_id) {
    bool valid;
    SeedEnvelope env = env_set_util::preprocess_env(peak_matrix, seed_env, mass_tol, &valid);
    if (!valid)
      return EnvSet();
    comp_peak_start_end_idx(peak_matrix, env.getPeakList(), mass_tol);
    std::vector<ExpEnvelope> exp_env_list;
    for (int idx = start_spec_id; idx < end_spec_id + 1; idx++) {
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix.get_row(idx), env, mass_tol);
      exp_env_list.push_back(exp_env);
    }
    EnvSet env_set = EnvSet(env, exp_env_list, start_spec_id, end_spec_id);
    return env_set;
  }
}
}