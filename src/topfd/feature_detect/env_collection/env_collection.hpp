//
// Created by abbash on 7/29/22.
//

#ifndef TOPPIC_ENV_COLLECTION_HPP
#define TOPPIC_ENV_COLLECTION_HPP

#include <memory>
#include <vector>
#include <string>
#include "topfd/feature_detect/env_set/env_set.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/spectrum/spectrum.hpp"
#include "topfd/feature_detect/util/utility_functions.hpp"
#include "common/xml/xml_dom_element.hpp"

namespace toppic {
class EnvCollection {
  public:
    EnvCollection();
    EnvCollection(const EnvCollection &ec);
    EnvCollection(const SeedEnvelope &env, const std::vector<EnvSet> &env_set_list, int min_charge, int max_charge, int start_spec_id, int end_spec_id);

    std::vector<std::vector<double>> get_seed_theo_map(PeakMatrix &peak_matrix, double snr);
    double comp_correlation();
    double comp_odd_even_log_ratio();
    void refine_mono_mass();
    double get_intensity(double snr, double noise_inte);
    double get_min_elution_time(spec_list &spectra_list);
    double get_max_elution_time(spec_list &spectra_list);
    double get_apex_elution_time(spec_list &spectra_list);
    double get_elution_length(spec_list &spectra_list);
    void remove_matrix_peaks(PeakMatrix &peak_matrix);
    void remove_peak_data(PeakMatrix &peak_matrix);
    std::vector<double> comp_exp_inte_sum_list();
    EnvSet get_seed_env_set();

    SeedEnvelope getSeedEnv() const { return seed_env_; }
    void setSeedEnv(const SeedEnvelope& seedEnv) { seed_env_ = seedEnv; }

    std::vector<EnvSet> getEnvSetList() const { return env_set_list_; }
    void setEnvSetList(const std::vector<EnvSet>& env_set_list) { for (auto & i : env_set_list) env_set_list_.push_back(i); }

    int getMinCharge() const { return min_charge_; }
    void setMinCharge(int minCharge) { min_charge_ = minCharge; }

    int getMaxCharge() const { return max_charge_; }
    void setMaxCharge(int maxCharge) { max_charge_ = maxCharge; }

    int getStartSpecId() const { return start_spec_id_; }
    void setStartSpecId(int startSpecId) { start_spec_id_ = startSpecId; }

    int getEndSpecId() const { return end_spec_id_; }
    void setEndSpecId(int endSpecId) { end_spec_id_ = endSpecId; }

    int getBaseSpecID() const { return seed_env_.getSpecId(); }
    double getMass() const { return seed_env_.getMass(); }
    std::vector<int> getChargeList() {
      std::vector<int> charge_list;
      for (auto es : env_set_list_)
        charge_list.push_back(es.getCharge());
      return charge_list;
    }

    std::vector<double> getExpInteSumList() const { return exp_inte_sum_list_; }
    void setExpInteSumList(const std::vector<double> &expInteSumList) {
      exp_inte_sum_list_.clear();
      for (auto & p : expInteSumList) exp_inte_sum_list_.push_back(p); }

    bool isEmpty(){
      if (seed_env_.isEmpty() && env_set_list_.empty() && min_charge_ == -1 && max_charge_ == -1 && start_spec_id_ == -1 && end_spec_id_ == -1)
        return true;
      return false;
    }

  private:
    SeedEnvelope seed_env_;
    std::vector<EnvSet> env_set_list_;
    int min_charge_;
    int max_charge_;
    int start_spec_id_;
    int end_spec_id_;
    std::vector<double> exp_inte_sum_list_;
};
}  // namespace toppic

#endif //TOPPIC_ENV_COLLECTION_HPP
