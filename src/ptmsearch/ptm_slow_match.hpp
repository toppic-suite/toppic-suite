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


#ifndef PROT_PTM_SLOW_MATCH_HPP_
#define PROT_PTM_SLOW_MATCH_HPP_

#include <memory>
#include <vector>

#include "base/proteoform.hpp"
#include "spec/prm_peak.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "oneptmsearch/diagonal_header.hpp"
#include "oneptmsearch/ptm_search_mng.hpp"
#include "oneptmsearch/basic_diag_pair.hpp"
#include "oneptmsearch/ps_align.hpp"

#include "ptmsearch/comp_shift_low_mem.hpp"

namespace prot {

class PtmSlowMatch {
 public:
  PtmSlowMatch(ProteoformPtr proteo_ptr,
               SpectrumSetPtr spectrum_set_ptr,
               AlignTypePtr align_type_ptr,
               CompShiftLowMem comp_shift,
               PtmSearchMngPtr mng_ptr);

  ProteoformPtr getProteoform(){return proteo_ptr_;};

  void init();

  PrsmPtr compute(AlignTypePtr align_type_ptr, int shift_num);

  void compute(AlignTypePtr type_ptr, PrsmPtrVec &prsm_ptrs);

 private:
  PtmSearchMngPtr mng_ptr_;
  ProteoformPtr proteo_ptr_;
  double prec_mono_mass_;
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  PrmMsPtrVec ms_six_ptr_vec_;
  ExtendMsPtrVec ms_three_ptr_vec_;
  AlignTypePtr align_type_ptr_;
  CompShiftLowMem comp_shift_;
  PSAlignPtr ps_align_ptr_;

  DiagonalHeaderPtrVec getNTermShiftListCommonHeaders();

  void addPrefixDiagonals(DiagonalHeaderPtrVec &common_header_ptrs,
                          DiagonalHeaderPtrVec &n_extend_header_ptrs);

  void addSuffixDiagonals(DiagonalHeaderPtrVec &common_header_ptrs,
                          DiagonalHeaderPtrVec &c_extend_header_ptrs);

  DiagonalHeaderPtrVec geneNTermShiftHeaders();
};

typedef std::shared_ptr<PtmSlowMatch> PtmSlowMatchPtr;
typedef std::vector<PtmSlowMatchPtr> PtmSlowMatchPtrVec;

} /* namespace prot */

#endif 
