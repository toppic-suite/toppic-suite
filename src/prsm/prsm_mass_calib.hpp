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

#ifndef TOPPIC_PRSM_PRSM_MASS_CALIB_HPP_
#define TOPPIC_PRSM_PRSM_MASS_CALIB_HPP_

#include <string>

#include "ms/calib/mass_calib_fit.hpp"
#include "para/prsm_para.hpp"

namespace toppic {

// The calibrants of toppic's built-in two-pass m/z recalibration (hidden
// --mass-calib; ms/calib, toppic_test experiment 152): the PrSMs of a
// zero-shift calibration search, read with their spectra, instead of the
// PrSM table that topcalib reads.
namespace prsm_mass_calib {

// The calibrants of <spectrum base name>.<input_file_ext>: the target PrSMs
// (protein names not starting with DECOY_) with an E-value <= max_evalue, no
// unexpected shift and no variable PTM. Precursors: the spectrum's precursor
// mass against the proteoform mass (mass_calib_fit::addPrecCalibrant).
// Fragments: the spectrum's masses within tol_ppm of the proteoform's
// theoretical fragment masses, as theo_peak_util generates them for the
// matching (the activation's ion types, the fixed and protein N-terminal
// modifications, the N-terminal label), the nearest one.
MassCalibrantSet readCalibrants(const PrsmParaPtr& prsm_para_ptr,
                                const std::string& input_file_ext,
                                double max_evalue, double tol_ppm);

// <base>_calib_ms2.msalign for <base>_ms2.msalign (the same directory): the
// corrected spectra that the main search reads
std::string calibSpFileName(const std::string& sp_file_name);

// <base>_mass_calib.json for <base>_ms2.msalign: the fitted model
std::string modelFileName(const std::string& sp_file_name);

}  // namespace prsm_mass_calib

}  // namespace toppic

#endif  // TOPPIC_PRSM_PRSM_MASS_CALIB_HPP_
