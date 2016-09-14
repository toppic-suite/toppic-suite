// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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
    }
    else if (type_ptr == RmBreakType::C_TERM) {
      if (!ion_ptr->getIonTypePtr()->isNTerm()) {
        cov = true;
      }
    }
    else if (type_ptr == RmBreakType::BOTH) {
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
  return cov_num/(double)total_num;
}

}
