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

#ifndef TOPPIC_MS_ENV_ENV_DETECT_HPP_
#define TOPPIC_MS_ENV_ENV_DETECT_HPP_

#include "ms/env/env.hpp"
#include "ms/env/env_para.hpp"
#include "ms/env/match_env.hpp"

namespace toppic {

namespace env_detect {

// min_inte: minimum intensity of an experimental peak counted as present in
// an envelope; min_ref_inte: minimum intensity of the reference (seed) peak;
// bound_inte: theoretical isotope peaks expected below this intensity (the
// noise level) are left out of the envelope.
MatchEnvPtr detectEnvByRefPeak(const PeakPtrVec& peak_list, int ref_peak,
                               int charge, double max_mass, double min_inte,
                               double min_ref_inte, double bound_inte,
                               const EnvParaPtr& env_para_ptr, bool is_full);

MatchEnvPtr detectEnvByMonoMass(const PeakPtrVec& peak_list, double mono_mass,
                                int charge, double min_inte, double bound_inte,
                                const EnvParaPtr& env_para_ptr, bool is_full);

MatchEnvPtr2D getCandidateEnv(const PeakPtrVec& peak_list, int max_charge,
                              double max_mass, double min_inte,
                              double min_ref_inte, double bound_inte,
                              const EnvParaPtr& env_para_ptr);
}  // namespace env_detect

}  // namespace toppic

#endif
