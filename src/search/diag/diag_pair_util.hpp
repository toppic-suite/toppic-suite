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

#ifndef TOPPIC_SEARCH_DIAG_DIAG_PAIR_UTIL_HPP_
#define TOPPIC_SEARCH_DIAG_DIAG_PAIR_UTIL_HPP_

#include "seq/proteoform.hpp"
#include "ms/spec/prm_peak.hpp"
#include "search/diag/diag_header.hpp"
#include "search/diag/diagonal.hpp"

namespace toppic {

namespace diag_pair_util {

DiagonalPtrVec geneDiagonalsWithEmptyList(const DiagHeaderPtrVec& header_ptr_vec,
                                          const PrmPeakPtrVec &prm_peaks,
                                          int group_spec_num, ProteoformPtr proteo_ptr);

DiagonalPtrVec geneDiagonalsWithoutEmptyList(const DiagHeaderPtrVec& header_ptr_vec,
                                             const PrmPeakPtrVec &prm_peaks,
                                             int group_spec_num, ProteoformPtr proteo_ptr);
}

} /* namespace toppic */

#endif 
