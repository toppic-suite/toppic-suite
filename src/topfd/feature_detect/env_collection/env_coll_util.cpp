//
// Created by abbash on 8/24/22.
//

#include "env_coll_util.hpp"
#include "topfd/feature_detect/env_set/env_set_util.hpp"

namespace toppic{
namespace env_coll_util{
  EnvCollection find_env_collection(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tole, int max_miss_env,
                                    int max_miss_charge, int max_miss_peak, int para_max_charge, double ratio_multi,
                                    double match_peak_tole) {
    SeedEnvelope env (seed_env);
    bool valid = env_set_util::check_valid_env_set(peak_matrix, env, mass_tole, match_peak_tole);
    if (!valid)
      return EnvCollection();
    int min_charge = 0, max_charge = 0;
    env_coll_util::get_charge_range(peak_matrix, env, mass_tole, max_miss_charge, para_max_charge, match_peak_tole, &min_charge, &max_charge);
    EnvSet top_peak_env_set = env_set_util::get_env_set_by_three_peaks(peak_matrix, env, mass_tole, max_miss_env);
    if (top_peak_env_set.isEmpty())
      return EnvCollection();
    int start_spec_id = top_peak_env_set.getStartSpecId();
    int end_spec_id = top_peak_env_set.getEndSpecId();
    double inte_ratio = top_peak_env_set.get_seed_inte_ratio() * ratio_multi;
    double base_inte = peak_matrix.get_min_inte();
    env.remove_low_inte_peaks(inte_ratio, base_inte);
    std::vector<EnvSet> env_set_list;
    for (int charge = min_charge; charge < max_charge + 1; charge++){
      SeedEnvelope cur_env = env.get_new_charge_env(charge);
      EnvSet env_set = env_set_util::find_env_set(peak_matrix, cur_env, mass_tole, start_spec_id, end_spec_id);
      if (!env_set.isEmpty()) {
        env_set.remove_all_non_consecutive_peaks(max_miss_peak);
        env_set_list.push_back(env_set);
      }
    }
    if (env_set_list.empty())
      return EnvCollection();
    EnvCollection env_coll = EnvCollection(env, env_set_list, min_charge, max_charge, start_spec_id, end_spec_id);
    return env_coll;
  }

  void get_charge_range(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tole, int max_miss_num,
                                  int para_max_charge, double match_peak_tole, int* return_min_charge, int* return_max_charge) {
    // search backward
    int min_charge = seed_env.getCharge();
    int charge = seed_env.getCharge() - 1;
    int miss_num = 0;
    while (charge >= 1) {
      SeedEnvelope env = seed_env.get_new_charge_env(charge);
      bool valid = env_set_util::check_valid_env_set(peak_matrix, env, mass_tole, match_peak_tole);
      if (!valid)
        miss_num = miss_num + 1;
      else {
        miss_num = 0;
        min_charge = charge;
      }
      if (miss_num >= max_miss_num)
        break;
      charge = charge - 1;
    }

    // search forward
    int max_charge = seed_env.getCharge();
    charge = seed_env.getCharge() + 1;
    while (charge <= para_max_charge) {
      SeedEnvelope env = seed_env.get_new_charge_env(charge);
      bool valid = env_set_util::check_valid_env_set(peak_matrix, env, mass_tole, match_peak_tole);
      if (!valid)
        miss_num = miss_num + 1;
      else {
        miss_num = 0;
        max_charge = charge;
      }
      if (miss_num >= max_miss_num)
        break;
      charge = charge + 1;
    }
    *return_min_charge = min_charge;
    *return_max_charge = max_charge;
  }

  SeedEnvelope select_best_seed_envelope(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tole, EnvBase env_base) {
    SeedEnvelope selected_envelope = SeedEnvelope();
    double max_corr = 0.5;
    for (int shift_num = -2; shift_num < 3; shift_num++) {
      SeedEnvelope s_env = seed_env.get_shifted_seed_envelope(env_base, shift_num);
      bool valid = env_set_util::check_valid_env_set(peak_matrix, s_env, mass_tole, 3);
      if (valid)
        env_set_util::comp_peak_start_end_idx(peak_matrix, s_env.getPeakList(), mass_tole);
      ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix.get_row(s_env.getSpecId()), s_env, mass_tole);
      std::vector<ExpPeak> exp_peak_list = exp_env.getExpEnvList();
      std::vector<double> exp_inte;
      for (auto p : exp_peak_list)
        if (!p.isEmpty())
          exp_inte.push_back(p.getInte());
        else
          exp_inte.push_back(0);
      double max_inte = *std::max(exp_inte.begin(), exp_inte.end());
      for (int i = 0; i < exp_inte.size(); i++)
           exp_inte[i] = exp_inte[i] / max_inte;

      std::vector<SimplePeak> theo_peak_list = s_env.getPeakList();
      std::vector<double> theo_inte;
      for (auto p : theo_peak_list)
        if (!p.isEmpty())
          theo_inte.push_back(p.getInte());
        else
          theo_inte.push_back(0);
      double corr = utility_functions::pearsonr(exp_inte, theo_inte);
      if (corr >= max_corr) {
        selected_envelope = s_env;
        max_corr = corr;
      }
    }
    return selected_envelope;
  }
}
}

