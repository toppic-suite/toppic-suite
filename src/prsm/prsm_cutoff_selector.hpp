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

#ifndef TOPPIC_PRSM_PRSM_CUTOFF_SELECTOR_HPP_
#define TOPPIC_PRSM_PRSM_CUTOFF_SELECTOR_HPP_

#include <string>

#include "common/base/mod.hpp"

namespace toppic {

namespace prsm_cutoff_selector {

// Write the PrSMs of <base_name>.<input_file_ext> that pass the cutoff to
// <base_name>.<output_file_ext>; base_name is the spectrum file name without
// its extension. cutoff_type: EVALUE (E-value), FDR (spectrum-level FDR),
// FORMFDR (spectrum- and proteoform-level FDR), PROTFDR (protein-level FDR)
// or FRAG (number of matched fragments, kept as the score). The proteoform
// of every written PrSM gets its classification (Proteoform::
// getClassification with its protein hits in db_file_name and var_mod_ptrs,
// the variable modifications of the search), unless it has one.
void process(const std::string& db_file_name, const std::string& base_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext, const std::string& cutoff_type,
             double cutoff_value, const ModPtrVec& var_mod_ptrs);
}

} /* namespace toppic */

#endif
