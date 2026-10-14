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

#ifndef TOPPIC_PRSM_PRSM_DEAMID_SELECTOR_HPP_
#define TOPPIC_PRSM_PRSM_DEAMID_SELECTOR_HPP_

#include <set>
#include <string>

namespace toppic {

// With deamidation as a variable PTM, the unmodified proteoform matched one
// isotope below the precursor mass (a zero-PTM PrSM at the -1 offset) and
// the deamidated proteoform matched at the precursor mass (a variable PTM
// PrSM) explain the same spectrum with masses 18.3 mDa apart. Their
// E-values come from PRM spectra of different precursor masses and do not
// decide between them (toppic_test experiment 117), so on each spectrum the
// two are compared by their matched fragments and precursor mass errors
// before the top PrSM selection (toppic_test experiment 123). toppic runs it
// only with the hidden --deamid-selector.
namespace prsm_deamid_selector {

// The abbreviations of the deamidations of a variable PTM file: its PTMs
// with Unimod id 7 or a monoisotopic mass within 0.001 Da of 0.984016 Da.
// Empty when the file has none.
std::set<std::string> getDeamidPtmNames(const std::string& var_ptm_file_name);

// Reads the PrSMs of input_file_ext (every PrSM of a spectrum together, in
// spectrum id order, as the E-value computation writes them) and writes
// them to output_file_ext without the losers below. On each spectrum, for
// every protein, start and end position that have both a zero-PTM PrSM (no
// unexpected shift, no variable PTM) and a deamidation PrSM (no unexpected
// shift; a variable PTM of deamid_ptm_names; the one with the best E-value
// if there are several), one of the two kinds is kept:
//   1. the one with more matched fragments when the numbers differ by 5 or
//      more;
//   2. otherwise the one whose precursor mass error (the precursor mass
//      moved by the whole isotopes nearest to the proteoform mass, minus
//      the proteoform mass) is closer to zero;
//   3. on an exact tie, the one with the better E-value.
// Every PrSM of the other kind in the group is removed. E-values are not
// changed.
void process(const std::string& spec_file_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext,
             const std::set<std::string>& deamid_ptm_names);

}  // namespace prsm_deamid_selector

}  // namespace toppic

#endif  // TOPPIC_PRSM_PRSM_DEAMID_SELECTOR_HPP_
