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

#ifndef TOPPIC_EXP_ENVELOPE_HPP
#define TOPPIC_EXP_ENVELOPE_HPP

#include <algorithm>
#include <vector>
#include "topfd/feature_detect/envelope/simple_peak.hpp"
#include "topfd/feature_detect/spectrum/exp_peak.hpp"

namespace toppic {
  class ExpEnvelope {
  public:
    ExpEnvelope();

    ExpEnvelope(int spec_id, std::vector<ExpPeak> peak_list);

    int get_match_peak_num(int base_idx);

    std::vector<double> get_inte_list();

    std::vector<double> get_pos_list();

    std::vector<double> get_non_empty_pos_list();

    int get_peak_num() { return peak_list_.size(); }

    ExpPeak get_peak(int idx) { return peak_list_[idx]; }

    int getSpecId() const { return spec_id_; }

    void setSpecId(int spec_id) { spec_id_ = spec_id; }

    std::vector<ExpPeak> getExpEnvList() { return peak_list_; }

    void setExpEnvList(const std::vector<ExpPeak> &peak_list);

    void setExpEnvListPeak(ExpPeak peak, int idx) { peak_list_[idx] = peak; }

    bool isEmpty();

  private:
    int spec_id_;
    std::vector<ExpPeak> peak_list_;
  };
}

#endif //TOPPIC_EXP_ENVELOPE_HPP
