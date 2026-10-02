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

#ifndef TOPPIC_COMMON_BASE_MOD_UTIL_HPP_
#define TOPPIC_COMMON_BASE_MOD_UTIL_HPP_

#include "common/base/mod.hpp"
#include "common/base/residue.hpp"

namespace toppic {

namespace mod_util {

ModPtrVec readModXml(const std::string& file_name);

std::vector<std::vector<std::string>> readModTxtForTsv(
    const std::string& file_name);

// The modifications of a modification file in four lists: [0] those that
// can occur at the N-terminus of a proteoform (Position N-term or any), [1]
// at its C-terminus (C-term or any), [2] anywhere (any) and [3] only at its
// N-terminus (N-term).
ModPtrVec2D readModTxt(const std::string& file_name);

ModPtrVec readAnywhereModTxt(const std::string& file_name);

// The modifications that can occur only at the N-terminus of a proteoform,
// on its first residue (Position N-term): list [3] of readModTxt.
ModPtrVec readNTermModTxt(const std::string& file_name);

std::vector<double> readModTxtToShiftList(const std::string& file_name);

ModPtrVec geneFixedModList(const std::string& str);

ResiduePtrVec geneResidueListWithMod(const ResiduePtrVec& residue_list,
                                     const ModPtrVec& fix_mod_list);

std::vector<double> getModMassVec(const ModPtrVec& var_mod_list);

}  // namespace mod_util

}  // namespace toppic

#endif
