// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#ifndef PROT_FEATURE_DECONV_DATA_HPP_
#define PROT_FEATURE_DECONV_DATA_HPP_

#include <memory>
#include <vector>

#include "spec/peak.hpp"

namespace prot {

class DeconvData {
 public:
  DeconvData(PeakPtrVec &peak_list, double max_mass, int max_charge, double win_size); 

	int getBgnPeak(int i) {return win_bgn_peaks_[i];}

	int getEndPeak(int i) {return win_end_peaks_[i];}

	int getIntervalPeakNum(int i) {return win_end_peaks_[i] - win_bgn_peaks_[i] + 1;}

	double getMaxMass() {return max_mass_;}

	int getMaxCharge() {return max_charge_;}

	PeakPtrVec getPeakList() {return peak_list_;}

	int getWinId(int i) {return win_ids_[i];}

	int getWinNum() {return win_num_;}

	void setMaxCharge(int max_chrg) {max_charge_ = max_chrg;}

	void setMaxMass(double mass) {max_mass_ = mass;}

 private:
	PeakPtrVec peak_list_;
  double max_mass_;
	int max_charge_;

	// the number of windows 
  int win_num_;
	// the length of each window 
  double win_size_;

	// the window id for each peak 
  std::vector<int> win_ids_;
	// the first peak of each window 
  std::vector<int> win_bgn_peaks_;
	// the last peak of each window 
  std::vector<int> win_end_peaks_;

  void initWinId();

  void initWinBgnEnd();
};

typedef std::shared_ptr<DeconvData> DeconvDataPtr;

}

#endif
