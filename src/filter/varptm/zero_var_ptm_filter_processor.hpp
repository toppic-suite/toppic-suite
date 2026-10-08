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

#ifndef TOPPIC_FILTER_VARPTM_ZERO_VAR_PTM_FILTER_PROCESSOR_HPP_
#define TOPPIC_FILTER_VARPTM_ZERO_VAR_PTM_FILTER_PROCESSOR_HPP_

#include "filter/mng/var_ptm_filter_mng.hpp"
#include "filter/mng/zero_ptm_filter_mng.hpp"

namespace toppic {

// The zero-shift filter and the variable PTM filter in one pass over the
// database blocks (toppic_test experiment 150). Both score the four
// zero-shift index files of a block (106 GB for the human database, more
// than the page cache holds), so the separate stages read them twice; here
// each block is loaded once by the zero-shift filter and the variable PTM
// filter scores on the same indexes. The candidates are those of the two
// stages run one after the other (the same toppic_zero_filter_* and
// toppic_var_filter_* files): each filter keeps its own spectrum sets,
// scores and candidate combiner.
namespace zero_var_ptm_filter_processor {

void process(const ZeroPtmFilterMngPtr& zero_mng_ptr,
             const VarPtmFilterMngPtr& var_mng_ptr);

}

} /* namespace toppic */

#endif
