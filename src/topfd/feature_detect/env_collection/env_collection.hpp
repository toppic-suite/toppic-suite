//
// Created by abbash on 7/29/22.
//
#include <memory>
#include <vector>
#include <string>

#include "topfd/feature_detect/env_set/env_set.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/spectrum/spectrum.hpp"
#include "topfd/feature_detect/util/utility_functions.hpp"

#include "common/xml/xml_dom_element.hpp"

#ifndef TOPPIC_ENV_COLLECTION_HPP
#define TOPPIC_ENV_COLLECTION_HPP

namespace toppic {
class EnvCollection {
    public:
      EnvCollection(SeedEnvelope env, std::vector<EnvSet> env_set_list, int min_charge, int max_charge, int start_spec_id, int end_spec_id){
        seed_env_ = env;
        for (auto & i : env_set_list) env_set_list_.push_back(i);
        min_charge_ = min_charge;
        max_charge_ = max_charge;
        start_spec_id_ = start_spec_id;
        end_spec_id_ = end_spec_id;
      }

    EnvCollection(){ seed_env_ = SeedEnvelope(); min_charge_ = -1; max_charge_ = -1; start_spec_id_ = -1; end_spec_id_ = -1; }

      double comp_correlation();
      double comp_odd_even_log_ratio();
      void refine_mono_mass();
      double get_intensity();
      double get_min_elution_time(spec_list spectra_list);
      double get_max_elution_time(spec_list spectra_list);
      double get_apex_elution_time(spec_list spectra_list);
      double get_elution_length(spec_list spectra_list);
      double remove_matrix_peaks(PeakMatrix peak_matrix);
      std::vector<double> comp_exp_inte_sum_list();


      SeedEnvelope getSeedEnv() const { return seed_env_; }
      void setSeedEnv(SeedEnvelope seedEnv) { seed_env_ = seedEnv; }

      std::vector<EnvSet> getEnvSetList() const { return env_set_list_; }
      void setEnvSetList(std::vector<EnvSet> env_set_list) { for (auto & i : env_set_list) env_set_list_.push_back(i); }

      int getMinCharge() const { return min_charge_; }
      void setMinCharge(int minCharge) { min_charge_ = minCharge; }

      int getMaxCharge() const { return max_charge_; }
      void setMaxCharge(int maxCharge) { max_charge_ = maxCharge; }

      int getStartSpecId() const { return start_spec_id_; }
      void setStartSpecId(int startSpecId) { start_spec_id_ = startSpecId; }

      int getEndSpecId() const { return end_spec_id_; }
      void setEndSpecId(int endSpecId) { end_spec_id_ = endSpecId; }

      std::vector<double> getExpInteSumList() const { return exp_inte_sum_list_; }
      void setExpInteSumList(std::vector<double> expInteSumList) { exp_inte_sum_list_ = expInteSumList; }

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
