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

#include <vector>
#include <algorithm>

#include "base/logger.hpp"
#include "base/base_algo.hpp"
#include "spec/theo_peak.hpp"
#include "spec/theo_peak_util.hpp"
#include "prsm/peak_ion_pair_util.hpp"

namespace prot {

namespace peak_ion_pair_util {

PeakIonPairPtrVec getMatchedPairs(const PeakIonPairPtrVec &pair_ptrs,
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

int getPeakIonPairNum(PeakIonPairPtrVec pairs) {
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

double computePairConverage(const PeakIonPairPtrVec &pair_ptrs, int begin,
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

PeakIonPairPtrVec findPairs(ExtendMsPtr ms_three_ptr,
                            TheoPeakPtrVec &theo_peak_ptrs,
                            int bgn, int end, double add_tolerance) {
  std::sort(theo_peak_ptrs.begin(), theo_peak_ptrs.end(), TheoPeak::cmpPosInc);
  std::vector<double> ms_masses = extend_ms::getExtendMassVec(ms_three_ptr);
  std::vector<double> theo_masses = theo_peak_util::getTheoMassVec(theo_peak_ptrs);

  PeakIonPairPtrVec pair_ptrs;
  size_t i = 0;
  size_t j = 0;
  while (i < ms_masses.size() && j < theo_masses.size()) {
    double deviation = ms_masses[i] - theo_masses[j];
    IonPtr ion_ptr = theo_peak_ptrs[j]->getIonPtr();
    double err = ms_three_ptr->getPeakPtr(i)->getOrigTolerance() + add_tolerance;
    if (ion_ptr->getPos() >= bgn && ion_ptr->getPos() <= end) {
      if (std::abs(deviation) <= err) {
        PeakIonPairPtr pair_ptr
            = std::make_shared<PeakIonPair>(ms_three_ptr->getMsHeaderPtr(),
                                            ms_three_ptr->getPeakPtr(i),
                                            theo_peak_ptrs[j]);
        pair_ptrs.push_back(pair_ptr);
      }
    }
    if (base_algo::increaseIJ(i, j, deviation, err, ms_masses, theo_masses)) {
      i++;
    } else {
      j++;
    }
  }
  return pair_ptrs;
}

/* parameter min_mass is necessary */
PeakIonPairPtrVec genePeakIonPairs(const ProteoformPtr &proteoform_ptr,
                                   const ExtendMsPtr &ms_three_ptr,
                                   double min_mass) {
  ActivationPtr activation_ptr
      = ms_three_ptr->getMsHeaderPtr()->getActivationPtr();

  TheoPeakPtrVec theo_peaks = theo_peak_util::geneProteoformTheoPeak(proteoform_ptr,
                                                                     activation_ptr,
                                                                     min_mass);

  return findPairs(ms_three_ptr, theo_peaks, 0, proteoform_ptr->getLen(), 0);
}

PeakIonPairPtrVec genePeakIonPairs(const ProteoformPtr &proteoform_ptr,
                                   const ExtendMsPtrVec &ms_ptr_vec,
                                   double min_mass) {
  PeakIonPairPtrVec pair_ptrs;
  for (size_t i = 0; i < ms_ptr_vec.size(); i++) {
    PeakIonPairPtrVec pair_ptr_tmp = genePeakIonPairs(proteoform_ptr, ms_ptr_vec[i],
                                                      min_mass);
    pair_ptrs.insert(pair_ptrs.end(), pair_ptr_tmp.begin(), pair_ptr_tmp.end());
  }
  return pair_ptrs;
}

double compMatchFragNum(const PeakIonPairPtrVec &pairs) {
  double match_fragment_num = 0;
  TheoPeakPtr prev_ion(nullptr);;
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i]->getTheoPeakPtr() != prev_ion) {
      prev_ion = pairs[i]->getTheoPeakPtr();
      match_fragment_num += pairs[i]->getRealPeakPtr()->getScore();
    }
  }
  return match_fragment_num;
}

double compMatchPeakNum(PeakIonPairPtrVec &pairs) {
  double match_peak_num = 0;
  std::sort(pairs.begin(), pairs.end(), PeakIonPair::cmpRealPeakPosInc);
  DeconvPeakPtr prev_deconv_peak(nullptr);
  LOG_DEBUG("total peak number " << pairs.size());
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i]->getRealPeakPtr()->getBasePeakPtr() != prev_deconv_peak) {
      prev_deconv_peak = pairs[i]->getRealPeakPtr()->getBasePeakPtr();
      match_peak_num += pairs[i]->getRealPeakPtr()->getScore();
    }
  }
  return match_peak_num;
}

}  // namespace peak_ion_pair_util

}  // namespace prot
