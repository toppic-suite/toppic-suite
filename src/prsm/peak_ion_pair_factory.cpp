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


#include <cmath>
#include <algorithm>
#include <vector>

#include "base/logger.hpp"
#include "base/base_algo.hpp"
#include "spec/extend_ms.hpp"
#include "spec/theo_peak.hpp"
#include "spec/theo_peak_factory.hpp"
#include "spec/theo_peak_util.hpp"
#include "prsm/peak_ion_pair_factory.hpp"

namespace prot {

PeakIonPairPtrVec PeakIonPairFactory::findPairs(ExtendMsPtr ms_three_ptr,
                                                TheoPeakPtrVec &theo_peak_ptrs,
                                                int bgn, int end, double add_tolerance) {
  std::sort(theo_peak_ptrs.begin(), theo_peak_ptrs.end(), TheoPeak::cmpPosInc);
  std::vector<double> ms_masses = ExtendMs::getExtendMassVec(ms_three_ptr);
  std::vector<double> theo_masses = theo_peak_util::getTheoMassVec(theo_peak_ptrs);

  PeakIonPairPtrVec pair_ptrs;
  size_t i = 0;
  size_t j = 0;
  while (i < ms_masses.size() && j < theo_masses.size()) {
    double deviation = ms_masses[i] - theo_masses[j];
    IonPtr ion_ptr = theo_peak_ptrs[j]->getIonPtr();
    double err = ms_three_ptr->getPeakPtr(i)->getOrigTolerance() 
        + add_tolerance;
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
PeakIonPairPtrVec PeakIonPairFactory::genePeakIonPairs(const ProteoformPtr &proteoform_ptr,
                                                       const ExtendMsPtr &ms_three_ptr,
                                                       double min_mass) {
  ActivationPtr activation_ptr
      = ms_three_ptr->getMsHeaderPtr()->getActivationPtr();

  TheoPeakPtrVec theo_peaks = TheoPeakFactory::geneProteoformTheoPeak(proteoform_ptr,
                                                                      activation_ptr,
                                                                      min_mass);

  return findPairs(ms_three_ptr, theo_peaks, 0, proteoform_ptr->getLen(), 0);
}

PeakIonPairPtrVec PeakIonPairFactory::genePeakIonPairs(const ProteoformPtr &proteoform_ptr,
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

}  // namespace prot
