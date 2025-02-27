//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_PTM_SEARCH_PTM_SLOW_MATCH_HPP_
#define TOPPIC_PTM_SEARCH_PTM_SLOW_MATCH_HPP_

#include <memory>
#include <vector>

#include "seq/proteoform.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "search/diag/diag_header.hpp"
#include "search/oneptmsearch/ptm_search_mng.hpp"
#include "search/oneptmsearch/ps_align.hpp"

namespace toppic {

class PtmSlowMatch {
 public:
  PtmSlowMatch(ProteoformPtr proteo_ptr,
               SpectrumSetPtr spectrum_set_ptr,
               ProteoformTypePtr align_type_ptr,
               PtmSearchMngPtr mng_ptr);

  ProteoformPtr getProteoform(){return proteo_ptr_;}

  void init();

  PrsmPtr compute(ProteoformTypePtr align_type_ptr, int shift_num);

  void compute(ProteoformTypePtr type_ptr, PrsmPtrVec &prsm_ptrs);

 private:
  PtmSearchMngPtr mng_ptr_;
  ProteoformPtr proteo_ptr_;
  double prec_mono_mass_;
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  PrmMsPtrVec ms_six_ptr_vec_;
  ExtendMsPtrVec ms_three_ptr_vec_;
  ProteoformTypePtr align_type_ptr_;
  PsAlignPtr ps_align_ptr_;

  DiagHeaderPtrVec getNTermShiftListCommonHeaders();

  void addPrefixDiagonals(DiagHeaderPtrVec &common_header_ptrs,
                          DiagHeaderPtrVec &n_extend_header_ptrs);

  void addSuffixDiagonals(DiagHeaderPtrVec &common_header_ptrs,
                          DiagHeaderPtrVec &c_extend_header_ptrs);

  DiagHeaderPtrVec geneNTermShiftHeaders();
};

typedef std::shared_ptr<PtmSlowMatch> PtmSlowMatchPtr;
typedef std::vector<PtmSlowMatchPtr> PtmSlowMatchPtrVec;

} /* namespace toppic */

#endif
