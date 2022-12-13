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

#ifndef TOPPIC_ENV_COLLECTION_HPP
#define TOPPIC_ENV_COLLECTION_HPP

#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <valarray>
#include "topfd/feature_detect/score/env_util.hpp"
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

    EnvCollection(const SeedEnvelope &env, const std::vector<EnvSet> &env_set_list, int min_charge, int max_charge,
                  int start_spec_id, int end_spec_id);

    std::vector<double> comp_exp_inte_sum_list();

    bool isEmpty();

    std::vector<int> getChargeList();

    void refine_mono_mass();

    double get_intensity(double snr, double noise_inte);

    double get_min_elution_time(spec_list &spectra_list) { return spectra_list[start_spec_id_].getRt(); }

    double get_max_elution_time(spec_list &spectra_list) { return spectra_list[end_spec_id_].getRt(); }

    double get_elution_length(spec_list &spectra_list);

    double get_apex_elution_time(spec_list &spectra_list);

    void remove_peak_data(PeakMatrix &peak_matrix);

    std::vector<std::vector<double>> get_seed_theo_map(PeakMatrix &peak_matrix, double snr);

    EnvSet get_seed_env_set();

    SeedEnvelope getSeedEnv() const { return seed_env_; }

    void setSeedEnv(const SeedEnvelope &seedEnv) { seed_env_ = seedEnv; }

    std::vector<EnvSet> getEnvSetList() const { return env_set_list_; }

    void setEnvSetList(const std::vector<EnvSet> &env_set_list);

    int getMinCharge() const { return min_charge_; }

    void setMinCharge(int minCharge) { min_charge_ = minCharge; }

    int getMaxCharge() const { return max_charge_; }

    void setMaxCharge(int maxCharge) { max_charge_ = maxCharge; }

    int getStartSpecId() const { return start_spec_id_; }

    void setStartSpecId(int startSpecId) { start_spec_id_ = startSpecId; }

    int getEndSpecId() const { return end_spec_id_; }

    void setEndSpecId(int endSpecId) { end_spec_id_ = endSpecId; }

    double getEcscore() const { return ecscore_; }

    void setEcscore(double ecscore) { ecscore_ = ecscore; }

    int getBaseSpecID() const { return seed_env_.getSpecId(); }

    double getMass() const { return seed_env_.getMass(); }

    std::vector<double> getExpInteSumList() const { return exp_inte_sum_list_; }

    void setExpInteSumList(const std::vector<double> &expInteSumList);

  private:
    SeedEnvelope seed_env_;
    std::vector<EnvSet> env_set_list_;
    int min_charge_;
    int max_charge_;
    int start_spec_id_;
    int end_spec_id_;
    double ecscore_ = -1;
    std::vector<double> exp_inte_sum_list_;
  };
}  // namespace toppic

#endif //TOPPIC_ENV_COLLECTION_HPP
