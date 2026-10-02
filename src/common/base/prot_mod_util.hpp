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

#ifndef TOPPIC_COMMON_BASE_PROT_MOD_UTIL_HPP_
#define TOPPIC_COMMON_BASE_PROT_MOD_UTIL_HPP_

#include <vector>

#include "common/base/prot_mod.hpp"
#include "common/base/residue.hpp"

namespace toppic {

namespace prot_mod_util {

bool allowMod(const ProtModPtr& prot_mod_ptr, const ResiduePtrVec& residues);

bool containMod(const ProtModPtrVec& prot_mod_ptr_vec,
                const ProtModPtr& prot_mod_ptr);

// The masses of the modifications that the protein N-terminal forms of
// prot_mod_ptr_vec put on the first residue of a proteoform starting at
// residue start of the protein with residues: the forms that start there,
// modify their first residue (acetylation) and apply to the protein. A
// proteoform N-terminal variable PTM of such a mass at such a start is
// that form.
std::vector<double> getFormNTermModMasses(const ProtModPtrVec& prot_mod_ptr_vec,
                                          const ResiduePtrVec& residues,
                                          int start);

}  // namespace prot_mod_util

}  // namespace toppic

#endif
