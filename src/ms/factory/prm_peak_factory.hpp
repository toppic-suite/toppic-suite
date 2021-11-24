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

#ifndef TOPPIC_MS_FACTORY_PRM_PEAK_FACTORY_HPP_
#define TOPPIC_MS_FACTORY_PRM_PEAK_FACTORY_HPP_

#include "para/peak_tolerance.hpp"
#include "ms/spec/prm_peak.hpp"

namespace toppic {

namespace prm_peak_factory {

PrmPeakPtr getZeroPeakPtr(int spec_id, double prec_mono_mass,
                          PeakTolerancePtr tole_ptr, double score);

PrmPeakPtr getPrecPeakPtr(int spec_id, double prec_mono_mass,
                          PeakTolerancePtr tole_ptr, double score);

}  // namespace prm_peak_factory

}  // namespace toppic

#endif
