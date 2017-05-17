// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_MASS_ONE_PTM_FILTER_H_
#define PROT_MASS_ONE_PTM_FILTER_H_

#include "base/proteoform.hpp"
#include "spec/prm_ms.hpp"
#include "prsm/simple_prsm.hpp"
#include "zeroptmfilter/mass_match.hpp"
#include "oneptmfilter/one_ptm_filter_mng.hpp"

namespace prot {

class MassOnePtmFilter {
 public:
  MassOnePtmFilter(const ProteoformPtrVec &proteo_ptrs,
                   OnePtmFilterMngPtr mng_ptr);
  void computeBestMatch(const PrmMsPtrVec &prm_ms_ptr_vec,
                        const PrmMsPtrVec &srm_ms_ptr_vec);

  SimplePrsmPtrVec getCompMatchPtrs() {return comp_match_ptrs_;}
  SimplePrsmPtrVec getPrefMatchPtrs() {return pref_match_ptrs_;}
  SimplePrsmPtrVec getSuffMatchPtrs() {return suff_match_ptrs_;}
  SimplePrsmPtrVec getInternalMatchPtrs() {return internal_match_ptrs_;}

 private:
  OnePtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;

  MassMatchPtr diag_index_ptr_;
  MassMatchPtr rev_diag_index_ptr_;
  MassMatchPtr term_index_ptr_;
  MassMatchPtr rev_term_index_ptr_;

  SimplePrsmPtrVec comp_match_ptrs_;
  SimplePrsmPtrVec pref_match_ptrs_;
  SimplePrsmPtrVec suff_match_ptrs_;
  SimplePrsmPtrVec internal_match_ptrs_;
};

typedef std::shared_ptr<MassOnePtmFilter> MassOnePtmFilterPtr;
} /* namespace prot */

#endif /* ONE_PTM_FILTER_H_ */
