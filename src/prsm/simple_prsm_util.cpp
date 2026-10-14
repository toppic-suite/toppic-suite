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

#include "prsm/simple_prsm_util.hpp"

#include <algorithm>

#include "prsm/simple_prsm.hpp"

namespace toppic {

namespace simple_prsm_util {

// Two candidates of a spectrum are the same simple PrSM when they have the
// same protein and the same sequence range (issue #348); the one with the
// best score is kept. Candidates without a range (the multiple-shift
// filter's) have the range -1 to -1, so for them the protein decides.
SimplePrsmPtrVec getUniqueMatches(SimplePrsmPtrVec match_ptrs) {
  std::sort(match_ptrs.begin(), match_ptrs.end(),
            SimplePrsm::cmpSeqRangeIncScoreDec);
  SimplePrsmPtrVec unique_match_ptrs;
  for (size_t i = 0; i < match_ptrs.size(); i++) {
    if (!unique_match_ptrs.empty()) {
      const SimplePrsmPtr& prev = unique_match_ptrs.back();
      const SimplePrsmPtr& cur = match_ptrs[i];
      if (cur->getSeqName() == prev->getSeqName() &&
          cur->getSpectrumId() == prev->getSpectrumId() &&
          cur->getSeqStartPos() == prev->getSeqStartPos() &&
          cur->getSeqEndPos() == prev->getSeqEndPos()) {
        continue;
      }
    }
    unique_match_ptrs.push_back(match_ptrs[i]);
  }

  return unique_match_ptrs;
}

}  // namespace simple_prsm_util

}  // namespace toppic
