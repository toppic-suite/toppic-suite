//
// Created by abbash on 8/24/22.
//

#include "env_coll_util.hpp"

namespace toppic{
namespace env_coll_util{
    EnvCollection find_env_collection(PeakMatrix& peak_matrix, SeedEnvelope& env, double mass_tole, int max_miss_env,
                                    int max_miss_charge, int max_miss_peak, int para_max_charge, double ratio_multi,
                                    double match_peak_tole, double snr, std::string file_name) {
    EnvSet top_peak_env_set = env_set_util::get_env_set(peak_matrix, env, mass_tole, max_miss_env);
      if (top_peak_env_set.isEmpty())
        return EnvCollection();
    top_peak_env_set.refine_feature_boundary();
    if (!env_set_util::check_valid_env_set_seed_env(peak_matrix, top_peak_env_set))
      return EnvCollection();

    double even_odd_peak_ratios = component_score::get_agg_odd_even_peak_ratio(top_peak_env_set);
    if (std::abs(even_odd_peak_ratios) > 0.4) {
      env = utility_functions::test_half_charge_state(peak_matrix, env, top_peak_env_set, even_odd_peak_ratios, mass_tole);
      if (env.isEmpty())
        return EnvCollection();
      EnvSet tmp_peak_env_set = env_set_util::get_env_set(peak_matrix, env, mass_tole, max_miss_env);
      if (!tmp_peak_env_set.isEmpty()) {
        top_peak_env_set = tmp_peak_env_set;
        top_peak_env_set.refine_feature_boundary();
        if (!env_set_util::check_valid_env_set_seed_env(peak_matrix, top_peak_env_set))
          return EnvCollection();
      }
      else
        return EnvCollection();
    }

    int start_spec_id = top_peak_env_set.getStartSpecId();
    int end_spec_id = top_peak_env_set.getEndSpecId();
    std::vector<EnvSet> env_set_list = get_charge_env_list(peak_matrix, env, top_peak_env_set, para_max_charge, mass_tole, max_miss_peak, max_miss_charge, ratio_multi, snr);
    env_set_list.push_back(top_peak_env_set);
    std::sort(env_set_list.begin(), env_set_list.end(), EnvSet::cmpCharge);

    int min_charge = env_set_list[0].getCharge();
    int max_charge = env_set_list[env_set_list.size()-1].getCharge();
    if (env_set_list.empty())
      return EnvCollection();

    /// write env_set
//    for (auto env_set : env_set_list)
//      write_out_files::write_env_set(peak_matrix, env_set, file_name);

    EnvCollection env_coll = EnvCollection(env, env_set_list, min_charge, max_charge, start_spec_id, end_spec_id);
    return env_coll;
  }

  std::vector<EnvSet> get_charge_env_list(PeakMatrix &peak_matrix, SeedEnvelope &env, EnvSet &top_peak_env_set, double para_max_charge,
                                          double mass_tole, int max_miss_peak, double max_miss_charge, double max_miss_env, double snr) {
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
//          env_set.remove_all_non_consecutive_peaks(max_miss_peak);
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
//          env_set.remove_all_non_consecutive_peaks(max_miss_peak);
          env_set_list.push_back(env_set);
        }
      }
      if (miss_num >= max_miss_charge)
        break;
    }
    std::sort(env_set_list.begin(), env_set_list.end(), EnvSet::cmpCharge);
    return env_set_list;
  }

  bool check_overlap(PeakMatrix &peak_matrix, EnvCollection &f, EnvCollection &env_coll, double time_tol) {
    spec_list spectra_list = peak_matrix.get_spectra_list();
    double start_rt = f.get_min_elution_time(spectra_list);
    double end_rt = f.get_max_elution_time(spectra_list);
    double feature_start_rt = env_coll.get_min_elution_time(spectra_list);
    double feature_end_rt = env_coll.get_max_elution_time(spectra_list);

    double start = -1;
    if (start_rt <= feature_start_rt)
      start = feature_start_rt;
    else
      start = start_rt;

    double end = -1;
    if (start > -1) {
      if (end_rt <= feature_end_rt)
        end = end_rt;
      else
        end = feature_end_rt;
    }

    bool status = false;
    if (end > - 1) {
      double overlapping_rt_range = end - start;
      if (overlapping_rt_range > 0) {
        double feature_rt_range = feature_end_rt - feature_start_rt;
        double feature_coverage = overlapping_rt_range / feature_rt_range;
        if (feature_coverage > time_tol)
          status = true;
      }
    }
    return status;
  }

  bool check_charge_state_distance(const std::vector<int> &parent_charge_states, int charge_state) {
    int min_charge_diff = 10000000;
    for (auto parent_charge_state : parent_charge_states) {
      int charge_diff = abs(charge_state - parent_charge_state);
      if (charge_diff < min_charge_diff)
        min_charge_diff = charge_diff;
    }
    if (min_charge_diff > 2)
      return false;
    return true;
  }

  void merge_overlapping_charge_features(EnvCollection &f, EnvCollection &env_coll) {
    std::vector<int> parent_charge_states = f.getChargeList();
    std::vector<EnvSet> esl = f.getEnvSetList();
    for (auto &feature : env_coll.getEnvSetList())
      if (std::count(parent_charge_states.begin(), parent_charge_states.end(), feature.getCharge()))
        esl.push_back(feature);
    f.setEnvSetList(esl);
  }

  bool check_in_existing_features(PeakMatrix &peakMatrix, EnvCollection &env_coll, std::vector<EnvCollection> &env_coll_list, double match_envelope_tolerance, double time_tol) {
    double mass_tol = match_envelope_tolerance * env_coll.getMass();
    std::vector<int> neighborFeatureList_charge_states = env_coll.getChargeList();
    std::vector<double> extended_masses = {env_coll.getMass() - 1.00235, env_coll.getMass(), env_coll.getMass() + 1.00235};
    int num_env_colls = env_coll_list.size();
    /// RT OVERLAPPPPP
    std::vector<int> selected_features;
    for (int i = 0; i < num_env_colls; i++) {
      if (check_overlap(peakMatrix, env_coll_list[i], env_coll, time_tol)) {
        double min_mass_diff = 100000000;
        for (auto ext_mass : extended_masses) {
          double mass_diff = std::abs(ext_mass - env_coll_list[i].getMass());
          if (mass_diff < min_mass_diff)
            min_mass_diff = mass_diff;
        }
        if (min_mass_diff < mass_tol)
          selected_features.push_back(i);
      }
    }

    bool status = true;
    bool overlap_charge = false;
    for (auto f_idx : selected_features) {
      EnvCollection f = env_coll_list[f_idx];
      std::vector<int> parent_charge_states = f.getChargeList();
      for (auto charge_state : neighborFeatureList_charge_states) {
        status = check_charge_state_distance(parent_charge_states, charge_state);
        if (std::count(parent_charge_states.begin(), parent_charge_states.end(), charge_state))
          overlap_charge = true;
      }

      if (status and !overlap_charge) {
        std::vector<EnvSet> esl = f.getEnvSetList();
        for (const auto& feature: env_coll.getEnvSetList()) {
          esl.push_back(feature);
          std::cout << "Envelope belongs to Feature: " << f_idx << std::endl;
          return true;
        }
        f.setEnvSetList(esl);
      }

      if (status and overlap_charge) {
        std::cout << "*********** Charge already exist: " << f_idx << std::endl;

        merge_overlapping_charge_features(f, env_coll);
        return true;
      }
    }
    return false;
  }

//  bool check_in_existing_features(PeakMatrix &peakMatrix, EnvCollection &env_coll, std::vector<EnvCollection> &env_coll_list, double match_envelope_tolerance, double time_tol) {
//    std::vector<int> neighborFeatureList_charge_states = env_coll.getChargeList();
//    double mass_tol = match_envelope_tolerance * env_coll.getMass();
//    std::vector<double> extended_masses;
//    for (int i = -2; i <= 2;  i++) {
//      double new_mass = env_coll.getMass() + (i * 1.00235);
//      extended_masses.push_back(new_mass);
//    }
//
//    /// RT OVERLAPPPPP
//    std::vector<int> selected_features;
//    for (int i = 0; i < env_coll_list.size(); i++) {
//      if (check_overlap(peakMatrix, env_coll_list[i], env_coll, time_tol)) {
//        double min_mass_diff = 100000000;
//        for (auto ext_mass : extended_masses) {
//          double mass_diff = std::abs(ext_mass - env_coll_list[i].getMass());
//          if (mass_diff < min_mass_diff)
//            min_mass_diff = mass_diff;
//        }
//        if (min_mass_diff < mass_tol)
//          selected_features.push_back(i);
//      }
//    }
//
//    bool status = true;
//    bool overlap_charge = false;
//    for (auto f_idx : selected_features) {
//      EnvCollection f = env_coll_list[f_idx];
//      std::vector<int> parent_charge_states = f.getChargeList();
//      for (auto charge_state : neighborFeatureList_charge_states) {
//        status = check_charge_state_distance(parent_charge_states, charge_state);
//        if (std::count(parent_charge_states.begin(), parent_charge_states.end(), charge_state))
//          overlap_charge = true;
//      }
//
//      if (status and !overlap_charge) {
//        std::vector<EnvSet> esl = f.getEnvSetList();
//        for (const auto& feature: env_coll.getEnvSetList()) {
//          esl.push_back(feature);
//          std::cout << "Envelope belongs to Feature: " << selected_features[0] << std::endl;
//          return true;
//        }
//        f.setEnvSetList(esl);
//      }
//
//      if (status and overlap_charge) {
//        std::cout << "*********** Charge already exist: " << selected_features[0] << std::endl;
//        merge_overlapping_charge_features(f, env_coll);
//        return true;
//      }
//    }
//    return false;
//  }
}
}