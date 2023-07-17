//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_TOPFD_ECSCORE_ECSCORE_PARA_HPP_
#define TOPPIC_TOPFD_ECSCORE_ECSCORE_PARA_HPP_

#include <vector>

#include "para/peak_tolerance.hpp"
#include "topfd/common/topfd_para.hpp"

namespace toppic {

class EcscorePara {
 public:
  EcscorePara(int frac_id, const std::string &file_name, 
              const std::string &resource_dir);

  EcscorePara(int frac_id, const std::string &file_name,
              const std::string &resource_dir, TopfdParaPtr para_ptr);

  std::vector<double> getExtendMasses(double mass);

  std::vector<double> getExtendOffsets() {return extend_offsets_;}

  std::vector<double> getSearchMasses(double mass);

  std::vector<double> getSearchOffsets() {return search_offsets_;}

  PeakTolerancePtr peak_tolerance_ptr_;

  // extend_offset is -2, -1, 0, 1, 2. We use it to search matched masses 
  // when a reference mass is given to find a feature. 
  std::vector<double> extend_offsets_;

  // search offset is  -3, -2, -1, 0, 1, 2, 3. We use it to search the precursor 
  // mass of an MS/MS spectrum to find a matched feature.
  std::vector<double> search_offsets_;

  double extend_min_mass_ = 5000;
  int intv_width_ = 500;
  int feature_num_ = 20000;
  int frac_id_;

  ///// params
  int para_max_charge_;
  int para_min_charge_;
  bool filter_neighboring_peaks_;
  double corr_tole_;

  double bin_size_ = 0.1;
  double neighbor_mass_tole_ = 0.01;
  double mass_tole_ = 0.008;
  int max_miss_env_ = 2;
  int max_miss_charge_ = 2;
  int max_miss_peak_ = 2;

  int min_match_peak_ = 2;
  int min_seed_match_peak_ = 2;
  int min_match_env_ = 2;

  double match_envelope_tolerance_ = 10E-6;
  double time_overlap_tole_ = 0.8;
  double even_odd_ratio_cutoff_ = 0.4;

  std::string file_name_;

};

typedef std::shared_ptr<EcscorePara> EcscoreParaPtr;

} /* namespace */

#endif 
