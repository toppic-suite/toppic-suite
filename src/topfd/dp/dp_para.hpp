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

#ifndef TOPPIC_TOPFD_DP_DP_PARA_HPP_
#define TOPPIC_TOPFD_DP_DP_PARA_HPP_

#include <memory>
#include <vector>

namespace toppic {
  class DpPara;
  typedef std::shared_ptr<DpPara> DpParaPtr;

  class DpPara {
  public:
    // DpPara() and DpPara(DpParaPtr) will be removed in version 1.7
    DpPara(){};
    DpPara(DpParaPtr dp_para_ptr){
      check_double_increase_ = dp_para_ptr->check_double_increase_;
      coexist_table_ = dp_para_ptr->coexist_table_;
      max_env_num_per_peak_ = dp_para_ptr->max_env_num_per_peak_;
      dp_env_num_ = dp_para_ptr->dp_env_num_;
      max_env_num_per_vertex_ = dp_para_ptr->max_env_num_per_vertex_;
      mz_tolerance_ = dp_para_ptr->mz_tolerance_;
      dp_window_size_ = dp_para_ptr->dp_window_size_;
      env_num_per_win_ = dp_para_ptr->env_num_per_win_;
    };

    DpPara(double mz_tolerance){mz_tolerance_ = mz_tolerance;};
  

    //dynamic programming window size
    double dp_window_size_ = 1.0;
    // Env assigned to 1 m/z intervals
    // number of envelopes per window 
    // use a small number of envelopes to speed up computation
    int env_num_per_win_ = 5;

    // Check double increasing when two envelopes overlap 
    bool check_double_increase_ = true;
    std::vector<std::vector<bool>> coexist_table_;

    // maximum number of envelopes sharing one peak 
    int max_env_num_per_peak_ = 2;
    // used in dpB to specify the number of output envelopes 
    int dp_env_num_ = 300;
    // maximum number of vertices per window 
    int max_env_num_per_vertex_ = 10;

    //mz tolerance used in computing msalign score
    //mz tolerance is fixed, NOT reduced for high charge state envelopes
    double mz_tolerance_ = 0.02;
  };
  
} /* namespace */

#endif 
