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


#ifndef PROT_PRSM_PEAK_ION_PAIR_FACTORY_HPP_
#define PROT_PRSM_PEAK_ION_PAIR_FACTORY_HPP_

#include "base/proteoform.hpp"
#include "spec/extend_ms.hpp"
#include "spec/theo_peak.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

class PeakIonPairFactory {
 public:

  static PeakIonPairPtrVec findPairs(ExtendMsPtr ms_three_ptr, TheoPeakPtrVec &theo_peak_ptrs, 
                                  int bgn, int end, double add_tolerance);

  static PeakIonPairPtrVec genePeakIonPairs(const ProteoformPtr &proteoform_ptr, 
                                            const ExtendMsPtr &ms_three_ptr, double min_mass);

  static PeakIonPairPtrVec genePeakIonPairs(const ProteoformPtr &proteoform_ptr, 
                                            const ExtendMsPtrVec &ms_ptr_vec, double min_mass);
};

}
#endif
