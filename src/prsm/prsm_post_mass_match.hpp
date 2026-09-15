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

#ifndef TOPPIC_PRSM_PRSM_POST_MASS_MATCH_HPP_
#define TOPPIC_PRSM_PRSM_POST_MASS_MATCH_HPP_

#include <string>

#include "para/prsm_para.hpp"

namespace toppic {

// Post mass matching (on by default in toppic; --disable-post-match turns it
// off): after the search, match the
// theoretical fragment masses of each PrSM that none of TopFD's deconvoluted
// masses matched against the centroided MS/MS peaks that TopFD stored in its
// SQLite database, using the method of MSPathFinderT (Informed-Proteomics):
// for every charge state from 1 to min(15, precursor charge + 2), the
// theoretical isotopic envelope of the fragment is placed in the spectrum and
// accepted if its most abundant isotopic peak is present, at least
// min_peak_num of its isotopic peaks are present, and the observed intensities
// have a Pearson correlation >= 0.7 or a Bhattacharyya distance <= 0.03 with
// the theoretical ones. Each accepted fragment adds one experimental
// monoisotopic mass to the spectrum, whose confidence score is the EnvCNN
// score of its envelope, computed as TopFD does for the deconvoluted masses.
namespace prsm_post_mass_match {

// Reads the PrSMs of <spectrum base name>.<input_file_ext> with their spectra,
// runs the post mass matching, and writes
//   <base>_post_ms2.msalign     the spectra with the added masses (the input
//                              msalign file is left untouched),
//   <base>_post_ms2.<output_file_ext>  the PrSMs with the recounted matched
//                              masses and fragments,
//   <base>.sqlite               the added masses as new rows of the ms2_env and
//                              ms2_env_peak tables,
//   <base>_post_ms2.feature    a copy of the TopFD feature file, so the
//                              downstream steps find it under the new name.
// Returns the name of the new spectrum file, which the downstream steps read
// instead of the input one. If the SQLite database is missing (TopFD run with
// -N), nothing is matched and the files are written unchanged. The EnvCNN
// model must have been loaded with onnx_env_cnn::initModel.
std::string process(const PrsmParaPtr& prsm_para_ptr,
                    const std::string& input_file_ext,
                    const std::string& output_file_ext, int min_peak_num);

}  // namespace prsm_post_mass_match

}  // namespace toppic

#endif
