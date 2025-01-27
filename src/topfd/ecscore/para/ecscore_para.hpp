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

#ifndef TOPPIC_TOPFD_ECSCORE_PARA_ECSCORE_PARA_HPP_
#define TOPPIC_TOPFD_ECSCORE_PARA_ECSCORE_PARA_HPP_

#include <vector>

#include "para/peak_tolerance.hpp"
#include "topfd/common/topfd_para.hpp"

namespace toppic {

class EcscorePara {
 public:
  EcscorePara(int frac_id, const std::string &file_name,
              int max_charge, int min_scan_num);

  double getPeakMzTole() {return peak_mz_tole_;}

  int getMinMatchPeakNumInTopThree() {return min_match_peak_;}

  int frac_id_;
  std::string file_name_;

  ///// params
  int para_min_charge_ = 1;
  int para_max_charge_ = 30;
  double seed_env_inte_corr_tole_cutoff_ = 0.5;

  double bin_size_ = 0.1;
  double neighbor_mz_tole_ = 0.01;
  double peak_mz_tole_ = 0.008;
  int max_miss_env_ = 2;
  int max_miss_charge_ = 2;

  // min matched peak in an envelope
  int min_match_peak_ = 2;
  // min number of scans for a feature
  int min_scan_num_ = 3;
  
  // used for match two feature masses 10 ppm 
  double match_feature_ppm_tolerance_ = 0.000010;
  double match_feature_time_overlap_tole_ = 0.8;
  double even_odd_ratio_cutoff_ = 0.4;


};

typedef std::shared_ptr<EcscorePara> EcscoreParaPtr;

} /* namespace */

#endif 
