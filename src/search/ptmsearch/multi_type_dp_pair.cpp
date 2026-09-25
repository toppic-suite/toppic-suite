// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "search/ptmsearch/multi_type_dp_pair.hpp"

#include <limits>

#include "search/oneptmsearch/path_type.hpp"

namespace toppic {

MultiTypeDpPair::MultiTypeDpPair(int x, int y, double pair_score, double diff,
                                 int order, int n_shift,
                                 const DiagHeaderPtr& header_ptr)
    : Pair(x, y) {
  diff_ = diff;
  pair_score_ = pair_score;
  order_ = order;
  header_ptr_ = header_ptr;
  prev_pair_ptrs_.resize(kProteoformTypeNum);
  scores_.resize(kProteoformTypeNum);
  path_types_.resize(kProteoformTypeNum);
  for (int t = 0; t < kProteoformTypeNum; t++) {
    for (int i = 0; i < n_shift + 1; i++) {
      prev_pair_ptrs_[t].push_back(nullptr);
      scores_[t].push_back(-std::numeric_limits<double>::max());
      path_types_[t].push_back(path_type::TYPE_NULL);
    }
  }
}

void MultiTypeDpPair::updateTable(int type_id, int s, double score,
                                  int path_type,
                                  const MultiTypeDpPairPtr& prev_pair_ptr) {
  scores_[type_id][s] = score;
  path_types_[type_id][s] = path_type;
  prev_pair_ptrs_[type_id][s] = prev_pair_ptr;
}

} /* namespace toppic */
