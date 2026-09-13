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

#ifndef TOPPIC_PRSM_PRSM_MATCH_NUM_RECOUNT_HPP_
#define TOPPIC_PRSM_PRSM_MATCH_NUM_RECOUNT_HPP_

#include <string>

#include "para/prsm_para.hpp"

namespace toppic {

namespace prsm_match_num_recount {

// The identification stages search the deconvoluted spectra after removing
// the masses below the EnvCNN score cutoff (SpPara::getEnvCnnCutoff()), so
// the matched mass and matched fragment numbers they store in the PrSMs
// count only the remaining masses. This step reads the PrSMs of
// <spectrum base name>.<input_file_ext> together with the full spectra
// (no cutoff), recomputes the two numbers against all masses, and writes the
// PrSMs to <base>.<output_file_ext>. The input must be sorted by spectrum id.
void process(const PrsmParaPtr& prsm_para_ptr,
             const std::string& input_file_ext,
             const std::string& output_file_ext);

}  // namespace prsm_match_num_recount

}  // namespace toppic

#endif
