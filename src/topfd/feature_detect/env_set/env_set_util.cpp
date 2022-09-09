#include <valarray>
#include <iostream>
#include "env_set_util.hpp"

namespace toppic {
namespace env_set_util {
    ExpPeak pick_exp_peak(PeakMatrix& peak_matrix, SimplePeak& seed_peak, int sp_id, double mass_tol) {
      // get peaks within mass tolerance
      ExpPeak result_peak = ExpPeak();
      double max_inte = -1000000;
      double pos = seed_peak.getPos();
//      std::cout << sp_id << ", " << pos << std::endl;
      //add check to see if the start and end are valid???? >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      for (int idx = seed_peak.getStartIdx(); idx < seed_peak.getEndIdx() + 1; idx++) {
        std::vector<ExpPeak> bin_peaks = peak_matrix.get_bin_peak(sp_id, idx);
        for (const auto& matrix_peak : bin_peaks) {
          double mass_diff = std::abs(pos - matrix_peak.getPos());
//          std::cout << "MD: " << mass_diff << ", " << matrix_peak.getPos() << ", " << bin_peaks.size() << std::endl;
          if ( mass_diff < mass_tol && matrix_peak.getInte() > max_inte) {
            result_peak = matrix_peak; ////////////////////// this is slow!!!!
            max_inte = matrix_peak.getInte();
          }
        }
      }
      return result_peak;
    }

    ExpEnvelope get_match_exp_env(PeakMatrix& peak_matrix, SeedEnvelope& seed_env, int sp_id, double mass_tol) {
      std::vector<ExpPeak> peak_list;
      std::vector<SimplePeak> peaks = seed_env.getPeakList();
      ExpPeak peak = ExpPeak();
      for (auto& seed_peak : peaks) {
        if (seed_peak.getStartIdx() > -1 && seed_peak.getEndIdx() > -1)
          peak = pick_exp_peak(peak_matrix, seed_peak, sp_id, mass_tol);
        peak_list.push_back(peak);
      }
      ExpEnvelope exp_env = ExpEnvelope(sp_id, peak_list);
      return exp_env;
    }

//  ExpPeak pick_exp_peak(PeakRow& row, SimplePeak& seed_peak, double mass_tol) {
//    // get peaks within mass tolerance
//    ExpPeak result_peak = ExpPeak();
//    double max_inte = -1000000;
//    double pos = seed_peak.getPos();
//    //add check to see if the start and end are valid???? >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//    std::vector<std::vector<ExpPeak>> rows = row.getRow();
//    for (int idx = seed_peak.getStartIdx(); idx < seed_peak.getEndIdx() + 1; idx++) {
//      for (const auto& matrix_peak : rows[idx]) {
//        double mass_diff = std::abs(pos - matrix_peak.getPos());
//        if ( mass_diff < mass_tol && matrix_peak.getInte() > max_inte) {
//          result_peak = matrix_peak; ////////////////////// this is slow!!!!
//          max_inte = matrix_peak.getInte();
//        }
//      }
//    }
//    return result_peak;
//  }
//
//  ExpEnvelope get_match_exp_env(PeakRow& peak_row, SeedEnvelope& seed_env, double mass_tol) {
//    std::vector<ExpPeak> peak_list;
//    std::vector<SimplePeak> peaks = seed_env.getPeakList();
//    ExpPeak peak = ExpPeak();
//    for (auto& seed_peak : peaks) {
//      if (seed_peak.getStartIdx() > -1 && seed_peak.getEndIdx() > -1)
//        peak = pick_exp_peak(peak_row, seed_peak, mass_tol);
//      peak_list.push_back(peak);
//    }
//    ExpEnvelope exp_env = ExpEnvelope(peak_row.getSpecID(), peak_list);
//    return exp_env;
//  }

  bool check_valid_env_set(PeakMatrix& peak_matrix, EnvSet& env_set) {
    bool valid = true;
    int start_spec_id = env_set.getStartSpecId();
    int end_spec_id = env_set.getEndSpecId();
    int start_idx = std::max(0, start_spec_id);
    int end_idx = std::min(end_spec_id, peak_matrix.get_spec_num() - 1);
    int elems = 0;
    std::vector<double> env_xic = env_set.getXicEnvIntes();
    for (int i = start_idx; i < end_idx; i++)
      if (env_xic[i] > 0) elems++;
    if (elems < 2) valid = false;
    return valid;
  }

   void remove_non_match_envs(std::vector<ExpEnvelope>& env_list, int refer_idx) {
    int idx = env_list.size() - 1;
    while (idx >= 0) {
      ExpEnvelope env = env_list[idx];
      if (env.get_match_peak_num(refer_idx) < 2) ///////////////////////////////////////////////////////////////////////////
        env_list.erase(env_list.begin() + idx);
      else
        return;
      idx = idx - 1;
    }
  }

  void comp_peak_start_end_idx(PeakMatrix& peak_matrix, SeedEnvelope &seed_env, double error_tole) {
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

  EnvSet get_env_set(PeakMatrix& peak_matrix, SeedEnvelope& env, double mass_tol, double max_miss_env) {
    // for checking if envelope should be extended.
    double snr = 3.0;
    double noise_inte_level = peak_matrix.get_min_inte();
    int base_idx = env.getSpecId();
    std::vector<double> theo_envelope_inte = env.get_inte_list();
    int refer_idx = std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) - theo_envelope_inte.begin();

    // search backward
    std::vector<ExpEnvelope> back_env_list;
    int idx = env.getSpecId();
    int miss_num = 0;
    while (idx >= 0) {
//      PeakRow row =  peak_matrix.get_row(idx);
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix, env, idx, mass_tol);
      /// process peaks based on inte
      std::vector<double> experimental_envelope_inte = exp_env.get_inte_list();
      double inte_ratio = env_utils::calcInteRatio_scan(theo_envelope_inte, experimental_envelope_inte);
      for (int i = 0; i < theo_envelope_inte.size(); i++) {
        double peak_inte = theo_envelope_inte[i];
        if ((inte_ratio * peak_inte) < (noise_inte_level * snr))
          exp_env.setExpEnvListPeak(ExpPeak(), i);
      }
      ///
      back_env_list.push_back(exp_env);
      if (exp_env.get_match_peak_num(refer_idx) < 2)
        miss_num = miss_num + 1;
      else
        miss_num = 0;
      if (miss_num >= max_miss_env)
        break;
      idx = idx - 1;
    }
    env_set_util::remove_non_match_envs(back_env_list, refer_idx);

    // search forward
    std::vector<ExpEnvelope> forw_env_list;
    idx = env.getSpecId() + 1;
    miss_num = 0;
    while (idx < peak_matrix.get_spec_num()) {
//      PeakRow row =  peak_matrix.get_row(idx);
//      ExpEnvelope exp_env = env_set_util::get_match_exp_env(row, env, mass_tol);
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix, env, idx, mass_tol);
      /// process peaks based on inte
      std::vector<double> experimental_envelope_inte = exp_env.get_inte_list();
      double inte_ratio = env_utils::calcInteRatio_scan(theo_envelope_inte, experimental_envelope_inte);
      for (int i = 0; i < theo_envelope_inte.size(); i++) {
        double peak_inte = theo_envelope_inte[i];
        if ((inte_ratio * peak_inte) < (noise_inte_level * snr))
          exp_env.setExpEnvListPeak(ExpPeak(), i);
      }
      ///
      forw_env_list.push_back(exp_env);
      if (exp_env.get_match_peak_num(refer_idx) < 2)
        miss_num = miss_num + 1;
      else
        miss_num = 0;
      if (miss_num >= max_miss_env)
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
    if ((end_spec_id - start_spec_id) < 2) return EnvSet();
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

  EnvSet find_env_set(PeakMatrix& peak_matrix, SeedEnvelope& env, double mass_tol, int start_spec_id, int end_spec_id) {
    double noise_inte_level = peak_matrix.get_min_inte();
    double snr = 3.0;
    ExpEnvelope empty_exp_env = ExpEnvelope();
    std::vector<double> theo_envelope_inte = env.get_inte_list();
    int refer_idx = std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) - theo_envelope_inte.begin();

    int base_idx = env.getSpecId();
//    int start_idx = std::max(base_idx - 1, 0);
//    int end_idx = std::min(base_idx + 1, peak_matrix.get_spec_num() - 1);

    int miss_num = 0;
    int max_miss_env = 2;

    std::vector<ExpEnvelope> back_env_list;
    for (int idx = base_idx; idx > start_spec_id; idx--) {
//      PeakRow row =  peak_matrix.get_row(idx);
//      ExpEnvelope exp_env = env_set_util::get_match_exp_env(row, env, mass_tol);
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix, env, idx, mass_tol);
      /// process peaks based on inte
      std::vector<double> experimental_envelope_inte = exp_env.get_inte_list();
      double inte_ratio = env_utils::calcInteRatio_scan(theo_envelope_inte, experimental_envelope_inte);
      for (int i = 0; i < theo_envelope_inte.size(); i++) {
        double peak_inte = theo_envelope_inte[i];
        if ((inte_ratio * peak_inte) < (noise_inte_level * snr))
          exp_env.setExpEnvListPeak(ExpPeak(), i);
      }
      back_env_list.push_back(exp_env);
      if (exp_env.get_match_peak_num(refer_idx) < 2)
        miss_num = miss_num + 1;
      else
        miss_num = 0;
      if (miss_num >= max_miss_env)
        break;
    }
    env_set_util::remove_non_match_envs(back_env_list, refer_idx);

    std::vector<ExpEnvelope> forw_env_list;
    for (int idx = base_idx; idx < end_spec_id; idx++) {
//      PeakRow row =  peak_matrix.get_row(idx);
//      ExpEnvelope exp_env = env_set_util::get_match_exp_env(row, env, mass_tol);
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix, env, idx, mass_tol);
      /// process peaks based on inte
      std::vector<double> experimental_envelope_inte = exp_env.get_inte_list();
      double inte_ratio = env_utils::calcInteRatio_scan(theo_envelope_inte, experimental_envelope_inte);
      for (int i = 0; i < theo_envelope_inte.size(); i++) {
        double peak_inte = theo_envelope_inte[i];
        if ((inte_ratio * peak_inte) < (noise_inte_level * snr))
          exp_env.setExpEnvListPeak(ExpPeak(), i);
      }
      forw_env_list.push_back(exp_env);
      if (exp_env.get_match_peak_num(refer_idx) < 2)
        miss_num = miss_num + 1;
      else
        miss_num = 0;
      if (miss_num >= max_miss_env)
        break;
    }
    env_set_util::remove_non_match_envs(forw_env_list, refer_idx);

    std::reverse(back_env_list.begin(), back_env_list.end());
    back_env_list.insert(back_env_list.end(), forw_env_list.begin(), forw_env_list.end());
    if (back_env_list.empty()) return EnvSet();

    start_spec_id = back_env_list[0].getSpecId();
    end_spec_id = back_env_list[back_env_list.size() - 1].getSpecId();
    if ((end_spec_id - start_spec_id) < 2) return EnvSet();

    EnvSet env_set = EnvSet(env, back_env_list, start_spec_id, end_spec_id);
    return env_set;
  }
//  EnvSet find_env_set(PeakMatrix& peak_matrix, SeedEnvelope& env, double mass_tol, int start_spec_id, int end_spec_id) {
//    double noise_inte_level = peak_matrix.get_min_inte();
//    double snr = 3.0;
//    ExpEnvelope empty_exp_env = ExpEnvelope();
//    std::vector<double> theo_envelope_inte = env.get_inte_list();
//    int refer_idx = std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) - theo_envelope_inte.begin();
//
//    int base_idx = env.getSpecId();
//    int start_idx = std::max(base_idx - 1, 0);
//    int end_idx = std::min(base_idx + 1, peak_matrix.get_spec_num() - 1);
////    std::cout << base_idx << " " << start_idx << " " << end_idx << std::endl;
//
//    std::vector<ExpEnvelope> exp_env_list;
//    for (int idx = start_spec_id; idx < end_spec_id + 1; idx++) {
//      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix.get_row(idx), env, mass_tol);
//
//      /// process peaks based on inte
//      std::vector<double> experimental_envelope_inte = exp_env.get_inte_list();
//      double inte_ratio = env_utils::calcInteRatio_scan(theo_envelope_inte, experimental_envelope_inte);
//      for (int i = 0; i < theo_envelope_inte.size(); i++) {
//        double peak_inte = theo_envelope_inte[i];
//        if ((inte_ratio * peak_inte) < (noise_inte_level * snr))
//          exp_env.setExpEnvListPeak(ExpPeak(), i);
//      }
//
//      if (exp_env.get_match_peak_num(refer_idx) > 2)
//        exp_env_list.push_back(exp_env);
//      else
//        exp_env_list.push_back(empty_exp_env);
//
//      /// check if valid
////      if ((start_idx <= idx) && (idx <= end_idx))
////        if (exp_env.get_match_peak_num(refer_idx) < 2)
////          return EnvSet();
//    }
//    if (exp_env_list.empty()) return EnvSet();
//    env_set_util::remove_non_match_envs(exp_env_list, refer_idx);
//    start_spec_id = exp_env_list[0].getSpecId();
//    end_spec_id = exp_env_list[exp_env_list.size() - 1].getSpecId();
//    if ((end_spec_id - start_spec_id) < 2) return EnvSet();
//
//    EnvSet env_set = EnvSet(env, exp_env_list, start_spec_id, end_spec_id);
//    return env_set;
//  }
}
}