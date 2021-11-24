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

#ifndef TOPPIC_STAT_TDGF_COMP_PROB_VALUE_ARRAY_HPP_
#define TOPPIC_STAT_TDGF_COMP_PROB_VALUE_ARRAY_HPP_

#include "common/base/residue_freq.hpp"
#include "ms/spec/prm_peak.hpp"
#include "prsm/prsm.hpp"
#include "para/peak_tolerance.hpp"
#include "stat/tdgf/comp_prob_value.hpp"

namespace toppic {

namespace comp_prob_value_array {

void compProbArray(CompProbValuePtr comp_prob_ptr, const ResFreqPtrVec &n_term_residue_ptrs, 
                   const PrmPeakPtrVec2D &peak_ptr_2d, const PrsmPtrVec &prsm_ptrs, bool strict,
                   double prob_prec_mass, PeakTolerancePtr tole_ptr, std::vector<double> &results);

}

}

#endif
