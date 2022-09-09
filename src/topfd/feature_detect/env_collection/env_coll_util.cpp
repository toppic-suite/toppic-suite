//
// Created by abbash on 8/24/22.
//

#include "env_coll_util.hpp"
#include "topfd/feature_detect/env_set/env_set_util.hpp"

namespace toppic{
namespace env_coll_util{
    EnvCollection find_env_collection(PeakMatrix& peak_matrix, SeedEnvelope& env, double mass_tole, int max_miss_env,
                                    int max_miss_charge, int max_miss_peak, int para_max_charge, double ratio_multi,
                                    double match_peak_tole, double snr) {
    EnvSet top_peak_env_set = env_set_util::get_env_set(peak_matrix, env, mass_tole, max_miss_env);
      if (top_peak_env_set.isEmpty())
        return EnvCollection();
    top_peak_env_set.refine_feature_boundary();
    if (!env_set_util::check_valid_env_set(peak_matrix, top_peak_env_set))
      return EnvCollection();

    int start_spec_id = top_peak_env_set.getStartSpecId();
    int end_spec_id = top_peak_env_set.getEndSpecId();
    std::vector<EnvSet> env_set_list = get_charge_env_list(peak_matrix, env, top_peak_env_set, para_max_charge, mass_tole, max_miss_peak, max_miss_charge, ratio_multi, snr);
    env_set_list.push_back(top_peak_env_set);
    std::sort(env_set_list.begin(), env_set_list.end(), EnvSet::cmpCharge);

    int min_charge = env_set_list[0].getCharge();
    int max_charge = env_set_list[env_set_list.size()-1].getCharge();
    if (env_set_list.empty())
      return EnvCollection();
    EnvCollection env_coll = EnvCollection(env, env_set_list, min_charge, max_charge, start_spec_id, end_spec_id);
    return env_coll;
  }

  std::vector<EnvSet> get_charge_env_list(PeakMatrix& peak_matrix, SeedEnvelope& env, EnvSet top_peak_env_set, double para_max_charge, double mass_tole, int max_miss_peak, double max_miss_charge, double max_miss_env, double snr) {
    int start_spec_id = top_peak_env_set.getStartSpecId();
    int end_spec_id = top_peak_env_set.getEndSpecId();

    std::vector<EnvSet> env_set_list;
    int charge = env.getCharge() - 1;
    int miss_num = 0;
    while (charge >= 1) {
      SeedEnvelope cur_env = env.get_new_charge_env(charge);
      env_set_util::comp_peak_start_end_idx(peak_matrix, cur_env, mass_tole);
      EnvSet env_set = env_set_util::find_env_set(peak_matrix, cur_env, mass_tole, start_spec_id, end_spec_id);
      charge = charge - 1;
      if (env_set.isEmpty()) {
        miss_num = miss_num + 1;
      }
      else {
        env_set.refine_feature_boundary();
        if (!env_set_util::check_valid_env_set(peak_matrix, top_peak_env_set))
          miss_num = miss_num + 1;
        else {
          miss_num = 0;
          env_set.remove_all_non_consecutive_peaks(max_miss_peak);
          env_set_list.push_back(env_set);
        }
      }
      if (miss_num >= max_miss_charge)
        break;
    }

    miss_num = 0;
    charge = env.getCharge() + 1;
    while (charge <= para_max_charge) {
      SeedEnvelope cur_env = env.get_new_charge_env(charge);
      env_set_util::comp_peak_start_end_idx(peak_matrix, cur_env, mass_tole);
      EnvSet env_set = env_set_util::find_env_set(peak_matrix, cur_env, mass_tole, start_spec_id, end_spec_id);
      charge = charge + 1;
      if (env_set.isEmpty()) {
        miss_num = miss_num + 1;
      }
      else {
        env_set.refine_feature_boundary();
        if (!env_set_util::check_valid_env_set(peak_matrix, env_set))
          miss_num = miss_num + 1;
        else {
          miss_num = 0;
          env_set.remove_all_non_consecutive_peaks(max_miss_peak);
          env_set_list.push_back(env_set);
        }
      }
      if (miss_num >= max_miss_charge)
        break;
    }
    std::sort(env_set_list.begin(), env_set_list.end(), EnvSet::cmpCharge);
    return env_set_list;
  }
}
}

