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

#include <vector>
#include <algorithm>

#include "base/logger.hpp"
#include "base/algorithm.hpp"
#include "spec/theo_peak.hpp"
#include "prsm/peak_ion_pair_util.hpp"

namespace prot {

PeakIonPairPtrVec PeakIonPairUtil::getMatchedPairs(const PeakIonPairPtrVec &pair_ptrs,
                                                   int spec_id, int peak_id) {
  PeakIonPairPtrVec selected_pair_ptrs;
  for (size_t i = 0; i < pair_ptrs.size(); i++) {
    if (pair_ptrs[i]->getMsHeaderPtr()->getId() == spec_id &&
        pair_ptrs[i]->getRealPeakPtr()->getBasePeakPtr()->getId() == peak_id) {
      selected_pair_ptrs.push_back(pair_ptrs[i]);
    }
  }
  return selected_pair_ptrs;
}

int PeakIonPairUtil::getPeakIonPairNum(PeakIonPairPtrVec pairs) {
  int match_peak_num = 0;
  DeconvPeakPtr prev_deconv_peak(nullptr);
  std::sort(pairs.begin(), pairs.end(), PeakIonPair::cmpRealPeakPosInc);
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i]->getRealPeakPtr()->getBasePeakPtr() != prev_deconv_peak) {
      prev_deconv_peak = pairs[i]->getRealPeakPtr()->getBasePeakPtr();
      match_peak_num += pairs[i]->getRealPeakPtr()->getScore();
    }
  }
  return match_peak_num;
}

double PeakIonPairUtil::computePairConverage(const PeakIonPairPtrVec &pair_ptrs, int begin,
                                             int end, RmBreakTypePtr type_ptr) {
  int total_num = end - begin  + 1;
  if (total_num <= 0) {
    return 0.0;
  }
  std::vector<bool> is_cov(total_num);
  for (size_t i  = 0; i < pair_ptrs.size(); i++) {
    IonPtr ion_ptr = pair_ptrs[i]->getTheoPeakPtr()->getIonPtr();
    bool cov = false;
    if (type_ptr == RmBreakType::N_TERM) {
      if (ion_ptr->getIonTypePtr()->isNTerm()) {
        cov = true;
      }
    } else if (type_ptr == RmBreakType::C_TERM) {
      if (!ion_ptr->getIonTypePtr()->isNTerm()) {
        cov = true;
      }
    } else if (type_ptr == RmBreakType::BOTH) {
      cov = true;
    }
    if (cov) {
      int pos = ion_ptr->getPos();
      if (pos >= begin && pos <= end) {
        is_cov[pos - begin] = true;
      }
    }
  }
  int cov_num = 0;
  for (size_t i = 0; i < is_cov.size(); i++) {
    if (is_cov[i]) {
      cov_num++;
    }
  }
  return cov_num / static_cast<double>(total_num);
}
}  // namespace prot
