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

#ifndef TOPPIC_STAT_COUNT_TDGF_UTIL_HPP_
#define TOPPIC_STAT_COUNT_TDGF_UTIL_HPP_

#include <vector>
#include <string>

#include "seq/proteoform.hpp"
#include "common/base/residue_freq.hpp"

namespace toppic {

namespace count_util {

void updateNTermResidueCounts(ResiduePtrVec &residue_list, std::vector<double> &counts,
                              const ProteoformPtrVec &mod_proteo_ptrs);

void updateResidueCounts(const ResiduePtrVec &residue_list,
                         std::vector<double> &counts,
                         ProteoformPtr prot_ptr);

ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list,
                              const std::vector<double> &counts);

int computeAvgLength(const ResFreqPtrVec &residue_ptrs, double convert_ratio);

}  // namespace tdgf_util

}  // namespace toppic

#endif
