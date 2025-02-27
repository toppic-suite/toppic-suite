//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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


#include <cmath>

#include "ms/spec/baseline_util.hpp" 
#include "ms/spec/peak_list_util.hpp" 
#include "topfd/deconv/deconv_data.hpp" 

namespace toppic {

DeconvData::DeconvData(PeakPtrVec peak_list, double max_mass, int max_charge,
                       double win_size): 
    peak_list_(peak_list),
    max_mass_(max_mass),
    max_charge_(max_charge),
    win_size_(win_size) {
      win_num_ = (int) std::ceil(peak_list_util::findMaxPos(peak_list_) / win_size_) + 2;
      initWinId();
      initWinBgnEnd();
    }

DeconvData::DeconvData(PeakPtrVec peak_list, double max_mass, int max_charge,
                       double win_size, bool estimate_min_inte, double sn_ratio): 
    peak_list_(peak_list),
    max_mass_(max_mass),
    max_charge_(max_charge),
    win_size_(win_size) {
      win_num_ = (int) std::ceil(peak_list_util::findMaxPos(peak_list_) / win_size_) + 2;
      initWinId();
      initWinBgnEnd();
      initMinInte(estimate_min_inte, sn_ratio);
    }

// initialize the window id 
void DeconvData::initWinId() {
  for (size_t i = 0; i < peak_list_.size(); i++) {
    win_ids_.push_back(std::floor(peak_list_[i]->getPosition() / win_size_));
  }
}

// initialize the first and last peak for each window 
void DeconvData::initWinBgnEnd() {
  int peak_num = peak_list_.size();
  win_bgn_peaks_.resize(win_num_, peak_num);
  win_end_peaks_.resize(win_num_, -1);
  win_bgn_peaks_[0] = 0;
  win_end_peaks_[win_num_ - 1] = peak_num - 1;
  for (int i = 0; i < peak_num; i++) {
    int id = win_ids_[i];
    if (i < win_bgn_peaks_[id]) {
      win_bgn_peaks_[id] = i;
    }
    if (i > win_end_peaks_[id]) {
      win_end_peaks_[id] = i;
    }
  }

  int bgn = peak_num;
  for (int i = win_num_ - 1; i >= 0; i--) {
    if (win_bgn_peaks_[i] == peak_num) {
      win_bgn_peaks_[i] = bgn;
    } else {
      bgn = win_bgn_peaks_[i];
    }
  }
  int end = -1;
  for (int i = 0; i < win_num_; i++) {
    if (win_end_peaks_[i] == -1) {
      win_end_peaks_[i] = end;
    } else {
      end = win_end_peaks_[i];
    }
  }
}

void DeconvData::initMinInte(bool estimate_min_inte, 
                             double sn_ratio) {
  min_inte_ = 0;
  min_ref_inte_ = 0;
  if (estimate_min_inte) {
    std::vector<double> intes;
    for (size_t i = 0; i < peak_list_.size(); i++) {
      intes.push_back(peak_list_[i]->getIntensity());
    }
    min_inte_ = baseline_util::getBaseLine(intes);
    min_ref_inte_ = min_inte_ * sn_ratio;
  }
}

}
