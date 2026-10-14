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
// table and SQL writers need them. Masses whose EnvCNN score is below
// env_cnn_cutoff are left out; the default 0 keeps every mass (the reader
// gives a mass without a score the score 1). A PrSM gets the spectra of its
// spectrum id; with second_sp_file_name (the residual spectra of toppic
// --second-precursor, whose first precursor is the second precursor of the
// scan) it gets, from either file, the spectra of its spectrum id whose
// first precursor has its precursor id. The PrSMs are returned in the
// file's order, without those whose spectra are not found.
PrsmPtrVec readPrsmsWithSpectra(const PrsmParaPtr& prsm_para_ptr,
                                const std::string& input_file_ext,
                                double env_cnn_cutoff = 0.0,
                                const std::string& second_sp_file_name = "");

}  // namespace prsm_reader_util

} /* namespace toppic */

#endif /* TOPPIC_PRSM_READER_HPP_ */
