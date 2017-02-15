// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <cmath>

#include "feature/raw_ms_util.hpp" 
#include "feature/deconv_data.hpp" 

namespace prot {

DeconvData::DeconvData(PeakPtrVec &peak_list, double max_mass, int max_charge,
                       double win_size): 
    peak_list_(peak_list),
    max_mass_(max_mass),
    max_charge_(max_charge),
    win_size_(win_size) {
      win_num_ = (int) std::ceil(RawMsUtil::findMaxPos(peak_list_) / win_size_) + 2;
      initWinId();
      initWinBgnEnd();
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

}
