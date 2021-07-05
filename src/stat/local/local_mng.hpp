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
           double max_ptm_mass,
           double min_ptm_mass,
           const std::string &input_file_ext,
           const std::string &output_file_ext);

  PrsmParaPtr prsm_para_ptr_;

  std::string input_file_ext_;
  std::string output_file_ext_;
  std::string residueModFileName_;

  double threshold_;

  double theta_, beta_;

  double min_mass_;

  double max_ptm_mass_;

  double min_ptm_mass_;

  double p1_, p2_;
  
  double desc_ratio_ = 0.67;

  int LEFT_SUP_LIMIT_ = 10;

  int RIGHT_SUP_LIMIT_ = 10;

  int DESC_MATCH_LIMIT_ = 5;
};

typedef std::shared_ptr<LocalMng> LocalMngPtr;

}  // namespace toppic

#endif /* TOPPIC_STAT_LOCAL_MNG_HPP_ */
