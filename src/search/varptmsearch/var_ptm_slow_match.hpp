//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_SEARCH_VAR_PTM_SEARCH_VAR_PTM_SLOW_MATCH_HPP_
#define TOPPIC_SEARCH_VAR_PTM_SEARCH_VAR_PTM_SLOW_MATCH_HPP_

#include <memory>
#include <vector>

#include "seq/proteoform.hpp"
#include "ms/spec/prm_peak.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "search/diag/diag_header.hpp"
#include "search/varptmsearch/var_ptm_search_mng.hpp"
#include "search/varptmsearch/var_ptm_align.hpp"

namespace toppic {

class VarPtmSlowMatch {
 public:
  VarPtmSlowMatch(ProteoformPtr proteo_ptr,
                  SpectrumSetPtr spectrum_set_ptr,
                  VarPtmSearchMngPtr mng_ptr);

  ProteoformPtr getProteoform(){return proteo_ptr_;};

  void init();

  bool isSuccessInit() {return success_init_;}

  PrsmPtr compute(); 

 private:
  VarPtmSearchMngPtr mng_ptr_;
  ProteoformPtr proteo_ptr_;
  double prec_mono_mass_;
  double prec_error_tole_;
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  PrmMsPtrVec ms_six_ptr_vec_;
  VarPtmAlignPtr var_ptm_align_ptr_;
  bool success_init_;

  DiagHeaderPtrVec geneVarPtmNTermShiftHeaders();
};

typedef std::shared_ptr<VarPtmSlowMatch> VarPtmSlowMatchPtr;
typedef std::vector<VarPtmSlowMatchPtr> VarPtmSlowMatchPtrVec;

} /* namespace toppic */

#endif 
