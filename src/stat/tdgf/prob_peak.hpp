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

#ifndef TOPPIC_STAT_TDGF_PROB_PEAK_HPP_
#define TOPPIC_STAT_TDGF_PROB_PEAK_HPP_

#include "ms/spec/base_peak_type.hpp"
#include "ms/spec/prm_peak.hpp"

namespace toppic {

// A PRM peak of the generating function: its integer mass and tolerance.
// strict: the strict tolerance of every peak (the separate mode, the PRM
// spectra at the PrSM's precursor). Otherwise the peak's N-strict-C-relax
// tolerance: the strict one (ppm x m) of an ORIGINAL peak (a b-ion mass m)
// and the relaxed one (ppm x (m + P), the precursor's error added) of a
// REVERSED peak (the complementary mass P - m of a y-ion, P the observed
// precursor) - unless strict_complement, when a REVERSED peak gets its strict
// tolerance too (PrsmPara::useEvalueStrictComplement, toppic hidden
// --evalue-complement-tolerance, issue #616).
class ProbPeak {
 public:
  ProbPeak(const PrmPeakPtr& peak_ptr, int spectrum_id, int height, bool strict,
           bool strict_complement, double convert_ratio);
  int mass_;
  int tolerance_;
  BasePeakTypePtr base_type_ptr_;
  int spectrum_id_;
  int mass_bgn_;
  int mass_end_;
  int table_bgn_;
  int table_end_;
};

}  // namespace toppic

#endif
