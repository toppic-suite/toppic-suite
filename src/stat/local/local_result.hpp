
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

#ifndef TOPPIC_STAT_LOCAL_LOCAL_SCORE_HPP_
#define TOPPIC_STAT_LOCAL_LOCAL_SCORE_HPP_

#include <vector>
#include <memory>

#include "common/base/ptm.hpp"
#include "seq/proteoform.hpp"

namespace toppic {

class LocalResult {
 public:
  int match_score_ = -1;
  std::vector<double> scr_vec_; 
  ProteoformPtr form_ptr_;
  PtmPtr ptm_ptr_;
};

typedef std::shared_ptr<LocalResult> LocalResultPtr;

}

#endif
