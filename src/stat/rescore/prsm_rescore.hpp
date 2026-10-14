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

#ifndef TOPPIC_STAT_RESCORE_PRSM_RESCORE_HPP_
#define TOPPIC_STAT_RESCORE_PRSM_RESCORE_HPP_

#include <limits>
#include <string>
#include <vector>

#include "para/prsm_para.hpp"
#include "prsm/prsm.hpp"

namespace toppic {

namespace prsm_rescore {

// The feature names, in the order of compFeatures; the retention time
// residual of the model (rescore_model) comes after them.
const std::vector<std::string>& getFeatureNames();

// The features of a PrSM with its full spectra attached
// (prsm_reader_util::readPrsmsWithSpectra): the PrSM's scores, shifts,
// N-terminal form and terminus, its spectrum and precursor feature, and its
// b / y (c / z) matches recomputed (peak_ion_pair_util::genePeakIonPairs) on
// the refined spectra of filtered_prsm_ptr, the same PrSM read with the
// EnvCNN cutoff (null: the match features are 0).
std::vector<double> compFeatures(const PrsmPtr& prsm_ptr,
                                 const PrsmPtr& filtered_prsm_ptr,
                                 double min_mass);

// The PrSMs with the maximum E-value (ExpectedValue::getMaxDouble(), 1e300:
// mostly PrSMs that the cluster protein correction gave their cluster's
// proteoform and that then matched too few fragments for an E-value). The
// rescoring cannot rank them: their E-value feature is clipped like any
// poor E-value, and the features they share with their cluster (feature
// intensity, precursor, ECScore) can give them high scores (toppic_test
// experiment 153: on a Bruker run 465 of them passed the rescored 1% FDR,
// 47 of them entrapment PrSMs). With enabled, they get the lowest rescore,
// so that the FDRs rank them last, as the E-value does.
struct EValueFloor {
  bool enabled = true;
  // false: they are left out of the training; true: the model is trained
  // on every PrSM, as without the floor, and only their scores are replaced
  bool train = false;
  // their rescore: the most negative double, below every score the model
  // gives
  double value = std::numeric_limits<double>::lowest();
};

// Rescores the PrSMs of <base>.<input_file_ext> (targets and decoys, one per
// spectrum; <base> the spectrum file's name without its extension) with
// rescore_model and writes them with their <rescore> to
// <base>.<output_file_ext>, in the input order. The PrSMs of each search
// (toppic, topmg) are rescored apart, as the FDR counts them apart; a
// search whose model cannot be trained gets the -log10 E-value as its
// score. The spectra are read with the EnvCNN cutoff of -F, the masses the
// search used; with second_sp_file_name (the residual spectra of toppic
// --second-precursor) a PrSM of a second precursor gets its residual
// spectrum from it (prsm_reader_util::readPrsmsWithSpectra). floor: the
// PrSMs with the maximum E-value (EValueFloor); a search without a model
// keeps their -log10 E-value, already the lowest.
void process(const PrsmParaPtr& prsm_para_ptr,
             const std::string& input_file_ext,
             const std::string& output_file_ext,
             const std::string& second_sp_file_name = "",
             const EValueFloor& floor = EValueFloor());

}  // namespace prsm_rescore

}  // namespace toppic

#endif
