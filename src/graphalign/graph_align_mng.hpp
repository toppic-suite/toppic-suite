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


#ifndef PROT_GRAPH_ALIGN_MNG_HPP_
#define PROT_GRAPH_ALIGN_MNG_HPP_

#include <string>

#include "prsm/prsm_para.hpp"

namespace prot {

class GraphAlignMng {
 public:
  GraphAlignMng(PrsmParaPtr prsm_para_ptr, 
                const std::string & var_mod_file_name,
                int n_unknown_shift, int max_known_mods,
                int proteo_graph_gap, double max_ptm_mass,
                int thread_num,
                const std::string & input_file_ext,
                const std::string & output_file_ext) {
    prsm_para_ptr_ = prsm_para_ptr;
    var_mod_file_name_ = var_mod_file_name;
    n_unknown_shift_ = n_unknown_shift;
    max_known_mods_ = max_known_mods;
    proteo_graph_gap_ = proteo_graph_gap;
    max_ptm_mass_ = max_ptm_mass;
    thread_num_ = thread_num;
    input_file_ext_ = input_file_ext;
    output_file_ext_ = output_file_ext;
  }

  PrsmParaPtr prsm_para_ptr_;

  std::string var_mod_file_name_;

  // set it to 1 for testing 
  double error_tolerance_ = 0.1;

  double max_ptm_sum_mass_ = 10000.00;

  double min_consistent_dist_ = 1.0;

  double convert_ratio_ = 274.335215;
  //double convert_ratio_ = 1.0;

  int proteo_graph_gap_;

  int getIntTolerance() {return std::ceil(error_tolerance_ * convert_ratio_);}
  int getIntMaxPtmSumMass() {return std::ceil(max_ptm_sum_mass_ * convert_ratio_);}
  int getIntMinConsistentDist() {return std::ceil(min_consistent_dist_ * convert_ratio_);};

  double align_prefix_suffix_shift_thresh_ = 300;

  double refine_prec_step_width_ = 0.005;

  int max_known_mods_ = 10;

  int n_unknown_shift_ =2;

  int prec_error_ = 0;

  int thread_num_ = 1;

  double max_ptm_mass_ = 500;

  std::string input_file_ext_;

  std::string output_file_ext_;
};

typedef std::shared_ptr<GraphAlignMng> GraphAlignMngPtr;

} /* namespace_prot */

#endif
