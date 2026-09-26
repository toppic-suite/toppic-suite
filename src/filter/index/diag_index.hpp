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

#ifndef TOPPIC_FILTER_INDEX_DIAG_INDEX_HPP_
#define TOPPIC_FILTER_INDEX_DIAG_INDEX_HPP_

#include "filter/mng/diag_filter_mng.hpp"
#include "filter/mng/one_ptm_filter_mng.hpp"

namespace toppic {

namespace diag_index {

// Writes the multi_ptm_index file of every block. When one_ptm_mng_ptr is
// given and its column limits (max_proteoform_mass_, filter_scale_) and index
// parameters equal mng_ptr's, the file the one-shift stage already wrote as
// one_ptm_diag_index is byte-identical to what this stage would build, and
// it is hard-linked (copied where links are not possible) instead of rebuilt.
void process(const DiagFilterMngPtr& mng_ptr,
             const OnePtmFilterMngPtr& one_ptm_mng_ptr = nullptr);

}

} /* namespace toppic */

#endif /* DIAG_INDEX_PROCESSOR_HPP_ */
