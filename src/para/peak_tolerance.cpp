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

#include "para/peak_tolerance.hpp"

#include <cmath>

namespace toppic {

PeakTolerance::PeakTolerance(double ppo) : ppo_(ppo) {}

double PeakTolerance::compRelaxErrorTole(double m1, double m2) const {
  return compStrictErrorTole(m1 + m2);
}

int PeakTolerance::getIntPpm() const {
  return static_cast<int>(std::round(ppo_ * 1000000));
}

double PeakTolerance::compStrictErrorTole(double mass) const {
  double tolerance = mass * ppo_;
  if (use_min_tolerance_ && tolerance < min_tolerance_) {
    tolerance = min_tolerance_;
  }
  return tolerance;
}

}  // namespace toppic
