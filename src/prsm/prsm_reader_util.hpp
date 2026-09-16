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

#ifndef TOPPIC_PRSM_PRSM_READER_UTIL_HPP_
#define TOPPIC_PRSM_PRSM_READER_UTIL_HPP_

#include <fstream>

#include "common/base/mod.hpp"
#include "para/prsm_para.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_str.hpp"
#include "seq/fasta_index_reader.hpp"

namespace toppic {

namespace prsm_reader_util {

PrsmStrPtrVec readAllPrsmStrs(const std::string& input_file_name);

PrsmStrPtrVec readAllPrsmStrsMatchSeq(const std::string& input_file_name);

PrsmPtrVec readAllPrsms(const std::string& prsm_file_name,
                        const FastaIndexReaderPtr& fasta_reader_ptr,
                        const ModPtrVec& fix_mod_list);

// Read the PrSMs of <spectrum base name>.<input_file_ext> and attach to each
// its deconvoluted spectra (from the msalign file) and refined spectra, as the
// table and SQL writers need them.
PrsmPtrVec readPrsmsWithSpectra(const PrsmParaPtr& prsm_para_ptr,
                                const std::string& input_file_ext);

}  // namespace prsm_reader_util

} /* namespace toppic */

#endif /* TOPPIC_PRSM_READER_HPP_ */
