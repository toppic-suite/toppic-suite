// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_FILTER_VARPTM_VAR_PTM_FILTER_PROCESSOR_HPP_
#define TOPPIC_FILTER_VARPTM_VAR_PTM_FILTER_PROCESSOR_HPP_

#include "filter/massmatch/filter_spec.hpp"
#include "filter/mng/var_ptm_filter_mng.hpp"
#include "filter/varptm/var_ptm_filter.hpp"
#include "prsm/simple_prsm_mem_merge.hpp"

namespace toppic {

namespace var_ptm_filter_processor {

// the spectrum sets of the stage (precursor offset 0, with the truncation
// N-terminal modification rounds), converted once
FilterSpecPtrVec readSpectra(const VarPtmFilterMngPtr& mng_ptr);

// scores the spectra against one block's filter and adds the candidates to
// the combiner (also used by zero_var_ptm_filter_processor)
void scoreBlock(const VarPtmFilterPtr& filter_ptr, int block_idx,
                const VarPtmFilterMngPtr& mng_ptr,
                const FilterSpecPtrVec& spectra, SimplePrsmMemMerge* combiner);

void process(const VarPtmFilterMngPtr& mng_ptr);

}

} /* namespace toppic */

#endif
