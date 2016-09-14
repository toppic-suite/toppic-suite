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


#include <cmath>

#include "base/logger.hpp"
#include "base/algorithm.hpp"
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
  std::vector<double> theo_masses = TheoPeakUtil::getTheoMassVec(theo_peak_ptrs);

  PeakIonPairPtrVec pair_ptrs;
  size_t i = 0;
  size_t j = 0;
  while (i < ms_masses.size() && j < theo_masses.size()) {
    double deviation = ms_masses[i] - theo_masses[j];
    IonPtr ion_ptr = theo_peak_ptrs[j]->getIonPtr();
    double err = ms_three_ptr->getPeakPtr(i)->getOrigTolerance() + add_tolerance;
    if (ion_ptr->getPos() >= bgn && ion_ptr->getPos() <= end) {
      if (std::abs(deviation) <= err) {
        PeakIonPairPtr pair_ptr = PeakIonPairPtr(new PeakIonPair(
                ms_three_ptr->getMsHeaderPtr(), ms_three_ptr->getPeakPtr(i), theo_peak_ptrs[j]));
        pair_ptrs.push_back(pair_ptr);
      }
    }
    if (increaseIJ(i, j, deviation, err, ms_masses, theo_masses)) {
      i++;
    } else {
      j++;
    }
  }
  return pair_ptrs;
}

/* parameter min_mass is necessary */
PeakIonPairPtrVec PeakIonPairFactory::genePeakIonPairs (const ProteoformPtr &proteoform_ptr, 
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

}
