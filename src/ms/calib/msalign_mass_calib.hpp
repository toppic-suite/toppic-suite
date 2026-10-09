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

#ifndef TOPPIC_MS_CALIB_MSALIGN_MASS_CALIB_HPP_
#define TOPPIC_MS_CALIB_MSALIGN_MASS_CALIB_HPP_

#include <string>

#include "ms/calib/mass_calib_model.hpp"

namespace toppic {

namespace msalign_mass_calib {

// Corrects the masses of a TopFD MS/MS msalign file with a calibration model
// (the shortcut of toppic_test experiment 152: the deconvoluted masses
// instead of the raw peaks, no TopFD rerun) and writes the corrected file:
// - precursors (PRECURSOR_MZ / PRECURSOR_MASS, ':'-separated lists with
//   PRECURSOR_CHARGE): the MS1 model at the precursor m/z and the spectrum's
//   RETENTION_TIME (minutes): m/z' = m/z (1 - e 1e-6),
//   mass' = mass - e 1e-6 (mass + z proton);
// - fragment masses: the MS/MS model at the fragment's m/z
//   ((mass + z proton) / z), the same mass formula.
// The corrected values are written with six decimals, as TopFD writes them;
// every other line is copied unchanged. Returns the number of spectra; exits
// with an error message when a file cannot be opened or a spectrum has no
// RETENTION_TIME.
int process(const MassCalibModel& model, const std::string& in_file_name,
            const std::string& out_file_name);

}  // namespace msalign_mass_calib

}  // namespace toppic

#endif  // TOPPIC_MS_CALIB_MSALIGN_MASS_CALIB_HPP_
