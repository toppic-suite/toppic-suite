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


#ifndef ZERO_PTM_FILTER_ZERO_PTM_FILTER_H_
#define ZERO_PTM_FILTER_ZERO_PTM_FILTER_H_

#include "base/proteoform.hpp"
#include "spec/extend_peak.hpp"
#include "spec/extend_ms.hpp"
#include "prsm/simple_prsm.hpp"
#include "zeroptmfilter/zero_ptm_filter_mng.hpp"
#include "zeroptmfilter/comp_shift.hpp"

namespace prot {

class ZeroPtmFilter {
 public:
  ZeroPtmFilter(const ProteoformPtrVec &proteo_ptrs,
                ZeroPtmFilterMngPtr mng_ptr);

  void computeBestMatch(const ExtendMsPtrVec &ms_ptr_vec);

  SimplePrsmPtrVec getCompMatchPtrs() {return comp_match_ptrs_;}

  SimplePrsmPtrVec getPrefMatchPtrs() {return pref_match_ptrs_;}

  SimplePrsmPtrVec getSuffMatchPtrs() {return suff_match_ptrs_;}

  SimplePrsmPtrVec getInternalMatchPtrs() {return internal_match_ptrs_;}

 private:
  ZeroPtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  CompShiftPtr index_ptr_;

  SimplePrsmPtrVec comp_match_ptrs_;
  SimplePrsmPtrVec pref_match_ptrs_;
  SimplePrsmPtrVec suff_match_ptrs_;
  SimplePrsmPtrVec internal_match_ptrs_;
};

typedef std::shared_ptr<ZeroPtmFilter> ZeroPtmFilterPtr;
} /* namespace prot */

#endif /* ZERO_PTM_FILTER_H_ */
