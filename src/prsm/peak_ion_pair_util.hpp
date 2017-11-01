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


#ifndef PROT_PRSM_PEAK_ION_PAIR_UTIL_HPP_
#define PROT_PRSM_PEAK_ION_PAIR_UTIL_HPP_

#include "spec/rm_break_type.hpp"
#include "spec/extend_peak.hpp"
#include "spec/theo_peak.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

class PeakIonPairUtil {
 public:
  static PeakIonPairPtrVec getMatchedPairs(const PeakIonPairPtrVec &pair_ptrs, 
                                           int spec_id, int peak_id);

  static int getPeakIonPairNum(PeakIonPairPtrVec pair_ptrs); 

  static double computePairConverage(const PeakIonPairPtrVec &pair_ptrs, int begin, 
                                     int end, RmBreakTypePtr type_ptr);

};

}
#endif
