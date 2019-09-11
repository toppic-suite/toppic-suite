//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "filter/massmatch/filter_protein.hpp"

namespace toppic {

FilterProtein::FilterProtein(int protein_id, int score):
    protein_id_(protein_id),
    score_(score) {
    }

inline bool cmpScore(const std::pair<int, int> &a, const std::pair<int, int> &b) {
  return a.second > b.second;
}

FilterProteinPtrVec FilterProtein::geneResults(std::vector<std::pair<int, int>> &single_type_results,
                                               int threshold, int single_type_num) {
  FilterProteinPtrVec prot_results;
  std::sort(single_type_results.begin(), single_type_results.end(), cmpScore);
  int output_num = 0;
  for (int i = 0; i < single_type_num; i++) {
    if (i >= static_cast<int>(single_type_results.size())) {
      break;
    }
    LOG_DEBUG("rank " << i  <<  " score " << single_type_results[i].second);
    if (single_type_results[i].second >= threshold) {
      output_num++;
    } else {
      break;
    }
  }
  for (int i = 0; i < output_num; i++) {
    int prot_id = single_type_results[i].first;
    int score = single_type_results[i].second;
    FilterProteinPtr prot_ptr = std::make_shared<FilterProtein>(prot_id, score);
    prot_results.push_back(prot_ptr);
  }
  return prot_results;
}

} /* namespace toppic */

