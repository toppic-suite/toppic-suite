//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_GRAPH_ALIGN_MNG_HPP_
#define PROT_GRAPH_ALIGN_MNG_HPP_

#include <string>

#include "prsm/prsm_para.hpp"

namespace prot {

class GraphAlignMng {
 public:
  GraphAlignMng(PrsmParaPtr prsm_para_ptr, 
                const std::string & var_mod_file_name,
                int n_unknown_shift,
                int max_known_mods,
                int proteo_graph_gap,
                int var_ptm_in_gap,
                double max_ptm_mass,
                int thread_num,
                const std::string & input_file_ext,
                const std::string & output_file_ext):
      prsm_para_ptr_(prsm_para_ptr),
      var_mod_file_name_(var_mod_file_name),
      n_unknown_shift_(n_unknown_shift),
      max_known_mods_(max_known_mods),
      proteo_graph_gap_(proteo_graph_gap),
      var_ptm_in_gap_(var_ptm_in_gap),
      max_ptm_mass_(max_ptm_mass),
      thread_num_(thread_num),
      input_file_ext_(input_file_ext),
      output_file_ext_(output_file_ext) {}

  PrsmParaPtr prsm_para_ptr_;

  std::string var_mod_file_name_;

  int n_unknown_shift_ = 0;

  int max_known_mods_ = 10;

  int proteo_graph_gap_ = 40;

  int var_ptm_in_gap_;

  // set it to 1 for testing 
  double error_tolerance_ = 0.1;

  double max_ptm_sum_mass_ = 10000.00;

  double min_consistent_dist_ = 1.0;

  double convert_ratio_ = 274.335215;

  int getIntTolerance() {return std::ceil(error_tolerance_ * convert_ratio_);}
  int getIntMaxPtmSumMass() {return std::ceil(max_ptm_sum_mass_ * convert_ratio_);}
  int getIntMinConsistentDist() {return std::ceil(min_consistent_dist_ * convert_ratio_);};

  double align_prefix_suffix_shift_thresh_ = 300;

  double refine_prec_step_width_ = 0.005;

  int prec_error_ = 0;

  double max_ptm_mass_ = 500;

  int thread_num_ = 1;

  std::string input_file_ext_;

  std::string output_file_ext_;
};

typedef std::shared_ptr<GraphAlignMng> GraphAlignMngPtr;

} /* namespace_prot */

#endif
