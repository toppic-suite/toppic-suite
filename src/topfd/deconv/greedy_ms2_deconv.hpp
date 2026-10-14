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

#ifndef TOPPIC_TOPFD_DECONV_GREEDY_MS2_DECONV_HPP_
#define TOPPIC_TOPFD_DECONV_GREEDY_MS2_DECONV_HPP_

#include <string>

#include "ms/env/match_env.hpp"
#include "ms/spec/peak.hpp"
#include "topfd/common/topfd_para.hpp"

namespace toppic {

// The greedy MS/MS deconvolution, --ms2-deconv greedy, TopFD's default
// (msdeconv, TopDIA's default, uses the candidate envelopes + dynamic
// programming of DeconvSingleSp for MS/MS spectra instead). It is the
// subtractive search of greedy_ms1_deconv adapted to fragment ions, in
// which the charge state of a peak is decided with all the charge states of
// its fragment mass:
//   1. take the most intense peak of the residual spectrum as the seed;
//   2. for every charge z up to the precursor's with a peak one isotope
//      spacing from the seed, align the averagine distribution's apex with
//      the seed, fit it to the residual intensities and extend the same mass
//      to the charges z-1, z-2, ... and z+1, z+2, ... while their envelopes
//      fit; the charge whose series explains the most intensity (minus the
//      theoretical intensity predicted where there is no peak, weighted by
//      the squared isotope cosine) wins, unless it looks like a harmonic of
//      charge z/2 or z/3;
//   3. choose the monoisotopic mass among +-2 isotope offsets by the cosine
//      with the isotope intensities summed over the series' charges;
//   4. with the hidden --ms2-max-shared-peak-num k (default -1, no
//      limit: the step is skipped), drop the envelopes of the mass that
//      share more than k peaks (the spectrum peaks both match) with an
//      envelope of the same charge already reported in the spectrum; a mass
//      with no envelope left is not reported;
//   5. subtract the theoretical envelopes of the mass at all its charges: a
//      peak much more intense than its theoretical intensity keeps the
//      difference, any other matched peak is used up;
//   6. repeat from 1 until no peak above the reference intensity is left.
// The parameters and their sources are in GreedyMs2Para (cpp).
namespace greedy_ms2_deconv {

// Deconvolute one MS/MS peak list (sorted by m/z) with fragment masses below
// max_mass and charges up to max_charge (deconvMsTwo's bounds from the
// precursors). The envelopes are matched to the peaks of peak_list, have
// their EnvCNN scores and are filtered and sorted as DeconvSingleSp's
// result. debug_name, when not empty, is the file the decisions of this
// spectrum are written to.
MatchEnvPtrVec deconv(const PeakPtrVec& peak_list,
                      const TopfdParaPtr& topfd_para_ptr, double max_mass,
                      int max_charge, const std::string& debug_name = "");

}  // namespace greedy_ms2_deconv

}  // namespace toppic

#endif
