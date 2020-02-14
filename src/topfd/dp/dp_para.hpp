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

#ifndef TOPPIC_TOPFD_DP_DP_PARA_HPP_
#define TOPPIC_TOPFD_DP_DP_PARA_HPP_

#include <memory>
#include <vector>

namespace toppic {

class DpPara {
 public:
  // DP algorithm
  // Check double increasing when two envelopes overlap 
  bool check_double_increase_ = true;
  std::vector<std::vector<bool>> coexist_table_;

  // maximum number of envelopes sharing one peak 
  int max_env_num_per_peak_ = 2;
  // used in dpB to specify the number of output envelopes 
  int dp_env_num_ = 300;
  // maximum number of vertices per window 
  int max_env_num_per_vertex_ = 10;
};

typedef std::shared_ptr<DpPara> DpParaPtr;

} /* namespace */

#endif 
