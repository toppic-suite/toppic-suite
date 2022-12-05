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

#ifndef TOPPIC_STAT_LOCAL_MNG_HPP_
#define TOPPIC_STAT_LOCAL_MNG_HPP_

#include <string>

#include "para/prsm_para.hpp"

namespace toppic {

class LocalMng {
 public:
  LocalMng(PrsmParaPtr prsm_para_ptr,
           double local_threshold,
           const std::string& residueModFileName,
           double min_ptm_mass,
           double max_ptm_mass,
           const std::string &input_file_ext,
           const std::string &output_file_ext);

  PrsmParaPtr prsm_para_ptr_;
  PeakTolerancePtr peak_tole_ptr_;

  std::string input_file_ext_;
  std::string output_file_ext_;
  std::string residueModFileName_;

  double threshold_;

  double min_mass_;

  double min_ptm_mass_;

  double max_ptm_mass_;

  // prob ratio between a matched fragment and an unmatched fragment
  double prob_ratio_ = 23.1434;

  // threshold for the deduction of matched fragments
  // when a PTM is localized
  double desc_ratio_ = 0.67;
  int DESC_MATCH_LIMIT_ = 5;

  // the range for checking possible starting position
  // of the N-terminus
  int n_term_range_ = 5;
  
  // the range for checking possible starting position
  // of the C-terminus
  int c_term_range_ = 5;

  //double theta_ = 0.994;  // the weight for known/unknown ptm

  //double beta_ = 0.8;     // the weight for one/two ptm

  //double p1_ = 0.915258; 

  //double p2_ = 21.1822; 
};

typedef std::shared_ptr<LocalMng> LocalMngPtr;

}  // namespace toppic

#endif /* TOPPIC_STAT_LOCAL_MNG_HPP_ */
