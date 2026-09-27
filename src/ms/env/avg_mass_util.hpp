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

#ifndef TOPPIC_MS_ENV_AVG_MASS_UTIL_HPP_
#define TOPPIC_MS_ENV_AVG_MASS_UTIL_HPP_

#include <string>

#include "ms/spec/deconv_ms.hpp"

namespace toppic {

// Conversions between the monoisotopic and the average masses of the
// deconvoluted masses in msalign files. The average mass of an ion is the
// intensity-weighted average mass of the theoretical isotopic distribution
// of its monoisotopic mass (EnvBase::convertMonoMassToAvgMass), so a
// conversion followed by its inverse gives the monoisotopic mass back.
// EnvBase::initBase must have been called.
namespace avg_mass_util {

// A copy of a deconvoluted MS/MS spectrum in which the precursor masses
// (and m/z values) and the fragment masses are average masses. The input
// spectrum and its header are left unchanged.
DeconvMsPtr getAvgMassMs(const DeconvMsPtr& ms_ptr);

// The name of the file convertToMonoMass writes for an average-mass
// msalign file: <base>_mono_ms2.msalign for <base>_ms2.msalign (for a
// file without the _ms2 suffix, <name>_mono.msalign).
std::string monoMassFileName(const std::string& avg_sp_file_name);

// Write the spectra of an average-mass msalign file (topfd --average-mass)
// with the precursor and fragment masses converted to monoisotopic masses
// to monoMassFileName(avg_sp_file_name), copy the TopFD feature file
// <base>.feature to the new name when it exists, and return the new file
// name. The parameter block of the input file is kept, with a comment line
// about the conversion added.
std::string convertToMonoMass(const std::string& avg_sp_file_name);

}  // namespace avg_mass_util

}  // namespace toppic

#endif
