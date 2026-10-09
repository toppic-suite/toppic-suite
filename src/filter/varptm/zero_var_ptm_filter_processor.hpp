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

#include "filter/mng/one_ptm_filter_mng.hpp"
#include "filter/mng/var_ptm_filter_mng.hpp"
#include "filter/mng/zero_ptm_filter_mng.hpp"

namespace toppic {

// The zero-shift, variable PTM and (with one_mng_ptr) one-shift filters as
// one pass over the database blocks: a block's four zero-shift index files
// are loaded once and scored by the three filters, whose candidates and
// output files are those of the separate stages (toppic_test experiment
// 150). The one-shift filter's files are the zero-shift ones cut at its
// column limit, so it scores the loaded indexes with the windows cut
// there, and takes the variable PTM filter's scores (the same mass lists)
// where they fit below the limit.
namespace zero_var_ptm_filter_processor {

void process(const ZeroPtmFilterMngPtr& zero_mng_ptr,
             const VarPtmFilterMngPtr& var_mng_ptr,
             const OnePtmFilterMngPtr& one_mng_ptr = nullptr);

}

} /* namespace toppic */

#endif
