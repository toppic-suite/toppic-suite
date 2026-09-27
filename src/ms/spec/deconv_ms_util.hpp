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

#ifndef TOPPIC_MS_SPEC_DECONV_MS_UTIL_HPP_
#define TOPPIC_MS_SPEC_DECONV_MS_UTIL_HPP_

#include <cstddef>

#include "ms/spec/deconv_ms.hpp"

namespace toppic {

namespace deconv_ms_util {

DeconvMsPtrVec getRefineMsPtrVec(const DeconvMsPtrVec& deconv_ms_ptr_vec,
                                 double new_prec_mass);

void keepTopPeaks(DeconvMsPtrVec& deconv_ms_ptr_vec, size_t peak_num);

// Combines the masses of one spectrum that agree within ppo, a RELATIVE
// tolerance (1e-5 for 10 ppm), into one mass. TopFD reports the same fragment
// mass once per charge state, so a spectrum carries the same mass several
// times; those copies are counted separately by the filters, the search and
// the E-value computation, which inflates the scores of targets and decoys
// alike.
//
// Each group is seeded by its most intense mass and takes every mass not yet
// used within ppo of the seed; the group is then represented by that seed,
// passed through unchanged. So every mass in the deduplicated list is a mass
// the file reports, with its own peak id, charge, score and intensity - which
// matters because peak ids are NOT renumbered: they are the line index of the
// mass in the msalign file and the key that joins a PrSM to an
// undeduplicated read of the same file at report time. Pointer identity is
// preserved throughout. Groups are emitted at the seed's position, so the
// list keeps the EnvCNN score order the file is written in. ppo <= 0 returns
// the input.
DeconvPeakPtrVec combineCloseMasses(const DeconvPeakPtrVec& peak_ptr_list,
                                    double ppo);

// The same, one spectrum at a time, never across spectra.
DeconvMsPtrVec combineCloseMasses(const DeconvMsPtrVec& deconv_ms_ptr_vec,
                                  double ppo);

}  // namespace deconv_ms_util

}  // namespace toppic

#endif
