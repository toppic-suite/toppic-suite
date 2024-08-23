// Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#include "topdia/pseudo_spec/pseudo_peak.hpp"

namespace toppic {
PseudoPeak::PseudoPeak(double mass, double monoMz, int charge,
                         double intensity, double score, double corr,
                         double shared_inte, int ms2_cycle_span,
                         double apexDiffScan, double rtLow, double rtHigh,
                         int ms2_apex_cycle) {
  mass_ = mass;
  mono_mz_ = monoMz;
  charge_ = charge;
  intensity_ = intensity;
  score_ = score;
  corr_ = corr;
  shared_inte_ = shared_inte;
  ms2_cycle_span_ = ms2_cycle_span;
  apex_diff_scan_ = apexDiffScan;
  rt_low_ = rtLow;
  rt_high_ = rtHigh;
  ms2_apex_cycle_ = ms2_apex_cycle;
}

PseudoPeak::PseudoPeak(const PseudoPeak &peaks) {
  mass_ = peaks.mass_;
  mono_mz_ = peaks.mono_mz_;
  charge_ = peaks.charge_;
  intensity_ = peaks.intensity_;
  score_ = peaks.score_;
  corr_ = peaks.corr_;
  shared_inte_ = peaks.shared_inte_;
  rank_ = peaks.rank_;
  ms2_cycle_span_ = peaks.ms2_cycle_span_;
  apex_diff_scan_ = peaks.apex_diff_scan_;
  rt_low_ = peaks.rt_low_;
  rt_high_ = peaks.rt_high_;
  ms2_apex_cycle_ = peaks.ms2_apex_cycle_;
  ms2_feature_idx_ = peaks.ms2_feature_idx_;
}
}  // namespace toppic
