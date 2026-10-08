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

#ifndef TOPPIC_MS_FACTORY_PRM_MS_UTIL_HPP_
#define TOPPIC_MS_FACTORY_PRM_MS_UTIL_HPP_

#include "ms/spec/prm_ms.hpp"
#include "para/peak_tolerance.hpp"
#include "para/sp_para.hpp"

namespace toppic {

namespace prm_ms_util {

// The filters' mass list: the integer masses round(x * scale) of the PRM
// peaks of the spectra with their tolerances, plus a zero and a precursor
// entry per spectrum, sorted by mass. mode decides which peaks of the same
// integer mass share an entry (SpPara::FilterMassListMode).
std::vector<std::pair<int, int> > getIntMassErrorList(
    const PrmMsPtrVec& prm_ms_ptr_vec, const PeakTolerancePtr& tole_ptr,
    double scale, bool n_strict, bool c_strict,
    SpPara::FilterMassListMode mode);

PrmPeakPtrVec getPrmPeakPtrs(const PrmMsPtrVec& prm_ms_ptr_vec,
                             const PeakTolerancePtr& tole_ptr);

}  // namespace prm_ms_util

}  // namespace toppic

#endif
