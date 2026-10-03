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

#ifndef TOPPIC_TOPFD_DECONV_GREEDY_MS1_DECONV_HPP_
#define TOPPIC_TOPFD_DECONV_GREEDY_MS1_DECONV_HPP_

#include <string>

#include "ms/env/match_env.hpp"
#include "ms/spec/peak.hpp"
#include "topfd/common/topfd_para.hpp"

namespace toppic {

// The greedy MS1 deconvolution (topfd --ms1-greedy-deconv), an alternative
// to the candidate envelopes + dynamic programming of DeconvSingleSp for MS1
// spectra. It is a THRASH-like subtractive search (Horn et al., JASMS 2000)
// in which the charge state of a peak is decided with all the charge states
// of its mass:
//   1. take the most intense peak of the residual spectrum;
//   2. for every charge z, align the averagine distribution's apex with the
//      peak, fit it to the residual intensities and extend the same mass to
//      the charges z-1, z-2, ... and z+1, z+2, ... while their envelopes fit;
//      the charge whose series explains the most intensity (minus the
//      theoretical intensity predicted where the spectrum has no peak) wins;
//   3. choose the monoisotopic mass among isotope offsets of that series;
//   4. subtract the theoretical envelopes of the mass at all its charges: a
//      peak much more intense than its theoretical intensity keeps the
//      difference, any other matched peak is used up;
//   5. repeat from 1 until no peak above the reference intensity is left.
namespace greedy_ms1_deconv {

// Deconvolute one MS1 peak list (sorted by m/z). The envelopes are matched
// to the peaks of peak_list, have their EnvCNN scores and are sorted as
// DeconvSingleSp's result. debug_name, when not empty, is the file the
// decisions of this spectrum are written to.
MatchEnvPtrVec deconv(const PeakPtrVec& peak_list,
                      const TopfdParaPtr& topfd_para_ptr,
                      const std::string& debug_name = "");

// Whether main's precursor-window step runs before the greedy search
// (research switch prec_win of TOPFD_GREEDY_PARA)
bool isPrecWinFirst();

}  // namespace greedy_ms1_deconv

}  // namespace toppic

#endif
