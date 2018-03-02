//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef PROT_SPEC_PRM_MS_HPP_
#define PROT_SPEC_PRM_MS_HPP_

#include <memory>
#include <vector>
#include <utility>

#include "spec/prm_peak.hpp"
#include "spec/ms.hpp"
#include "spec/peak_tolerance.hpp"

namespace prot {

typedef std::shared_ptr<Ms<PrmPeakPtr> > PrmMsPtr;

typedef std::vector<PrmMsPtr> PrmMsPtrVec;

namespace prm_ms {

std::vector<std::pair<int, int> > getIntMassErrorList(const PrmMsPtrVec &prm_ms_ptr_vec,
                                                      PeakTolerancePtr tole_ptr,
                                                      double scale, bool n_strict, bool c_strict);

PrmPeakPtrVec getPrmPeakPtrs(const PrmMsPtrVec &prm_ms_ptr_vec, PeakTolerancePtr tole_ptr);

}  // namespace prm_ms

}  // namespace prot

#endif
