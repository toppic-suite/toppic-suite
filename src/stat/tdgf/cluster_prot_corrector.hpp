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

#ifndef TOPPIC_STAT_TDGF_CLUSTER_PROT_CORRECTOR_HPP_
#define TOPPIC_STAT_TDGF_CLUSTER_PROT_CORRECTOR_HPP_

#include <string>

#include "para/prsm_para.hpp"
#include "stat/tdgf/tdgf_mng.hpp"

namespace toppic {

namespace cluster_prot_corrector {

// Reads the PrSMs of <spectrum base name>.<input_file_ext> (clustered into
// proteoforms, i.e. after prsm_feature_cluster / prsm_simple_cluster and
// prsm_prot_cluster). In every proteoform cluster, the PrSMs whose protein
// differs from the protein of the cluster's PrSM with the lowest E-value
// (ties: the lowest sequence name) are the same proteoform mass matched to
// another protein; each of them is given the best PrSM's proteoform (and
// its mass as adjusted precursor mass), its matched experimental and
// theoretical fragment numbers are recounted on its spectra read with the
// search's EnvCNN score cutoff (SpPara::getEnvCnnCutoff, so that the numbers
// are comparable with the other PrSMs' ones), and its E-value is recomputed
// with the generating-function method of tdgf_mng_ptr on those spectra. All
// PrSMs are written, in the input order, to <base>.<output_file_ext>.
// Returns the number of corrected PrSMs. The fragment deduplication
// tolerance (SpPara::getFragDedupPpo, 0 = off) combines close fragment
// masses for the E-value computation only: the matched fragment numbers are
// recounted here as well, and those must come from the masses the file
// reports.
int process(const PrsmParaPtr& prsm_para_ptr, const TdgfMngPtr& tdgf_mng_ptr,
            const std::string& input_file_ext,
            const std::string& output_file_ext);

}  // namespace cluster_prot_corrector

}  // namespace toppic

#endif
