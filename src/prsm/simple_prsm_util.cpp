//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include <algorithm>

#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_util.hpp"

namespace prot {

namespace simple_prsm_util {

SimplePrsmPtrVec getUniqueMatches(SimplePrsmPtrVec match_ptrs) {
  std::sort(match_ptrs.begin(), match_ptrs.end(), SimplePrsm::cmpNameIncScoreDec);
  SimplePrsmPtrVec unique_match_ptrs;
  std::string prev_name = "";
  int prev_id = -1;
  for(size_t i = 0; i < match_ptrs.size(); i++){
    std::string cur_name = match_ptrs[i]->getSeqName();
    int cur_id = match_ptrs[i]->getSpectrumId();
    if (cur_name != prev_name || cur_id != prev_id) {
      unique_match_ptrs.push_back(match_ptrs[i]);
      prev_name = cur_name;
      prev_id = cur_id;
    }
  }

  return unique_match_ptrs;
}

} // namespace simple_prsm_util

} // namespace prot 
