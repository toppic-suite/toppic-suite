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

#ifndef TOPPIC_PRSM_PRSM_SECOND_PREC_HPP_
#define TOPPIC_PRSM_PRSM_SECOND_PREC_HPP_

#include <string>

#include "para/prsm_para.hpp"

namespace toppic {

// The second precursor search of toppic --second-precursor: TopFD lists
// every MS1 feature of an MS/MS spectrum's isolation window as a precursor,
// and the search uses the first. After the first pass the spectra with two
// or more precursors are searched again with the second precursor, without
// the masses the first pass's accepted PrSM explains (TopMPI's search of a
// multiplexed spectrum: Wang, Xiong and Liu, Proteomics 2025), and the two
// passes' PrSMs are reported together, a spectrum with up to one PrSM per
// precursor.
namespace prsm_second_prec {

// Writes the residual spectra: of the spectra of ori_sp_file_name (the
// TopFD msalign file the first pass searched) with two or more precursors,
// the precursor list with the second precursor first (ids renumbered), the
// masses without those that the spectrum's accepted PrSM of the first pass
// matches (record_file_name: the first pass's record PrSMs; accepted:
// a target with fdr <= cutoff_value when cutoff_type is FDR, or with an
// E-value <= cutoff_value otherwise; matched: within the fragment tolerance
// of one of the proteoform's fragment masses, as it is or one isotope
// lower or higher). Returns the number of spectra written.
int writeResidualSpectra(const PrsmParaPtr& prsm_para_ptr,
                         const std::string& record_file_name,
                         const std::string& ori_sp_file_name,
                         const std::string& cutoff_type, double cutoff_value,
                         const std::string& resid_sp_file_name);

// Writes the spectra of first_sp_file_name followed by those of
// second_sp_file_name with their spectrum ids (and MS1 ids) raised by
// SpPara::getMaxSpecNumPerFile(); the feature ids are kept.
void mergeSpectra(const std::string& first_sp_file_name,
                  const std::string& second_sp_file_name,
                  const std::string& out_sp_file_name);

// Merges the record PrSMs of the two passes into one set, one unit per
// spectrum and precursor: a second-pass PrSM whose protein (decoy prefix
// aside) or precursor feature a first-pass PrSM of the spectrum already has
// is kept only when its E-value is better, in place of that PrSM. The
// second-pass PrSMs get the spectrum ids of mergeSpectra, precursor id 1
// and the merged spectrum file's name (sp_short_name); the first pass's
// keep theirs. The output is in spectrum id order.
void mergePrsms(const std::string& first_prsm_file_name,
                const std::string& second_prsm_file_name,
                const std::string& out_prsm_file_name,
                const std::string& sp_short_name);

}  // namespace prsm_second_prec

}  // namespace toppic

#endif
