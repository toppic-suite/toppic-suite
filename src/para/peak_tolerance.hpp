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

#ifndef TOPPIC_PARA_PEAK_TOLERANCE_HPP_
#define TOPPIC_PARA_PEAK_TOLERANCE_HPP_

#include <memory>
#include <vector>

namespace toppic {

class PeakTolerance {
 public:
  explicit PeakTolerance(double ppo);

  // The error tolerance of a mass: mass x ppo_, or the Dalton tolerance when
  // one is set; at least min_tolerance_ when use_min_tolerance_. Mass 0 (the
  // zero PRM peak's N-terminal mass) carries no measurement error: it keeps
  // the minimum tolerance with a Dalton tolerance too, as under any ppo_.
  double compStrictErrorTole(double mass) const;

  // The error tolerance of a precursor mass in the precursor mass tests of
  // the filters, the searches and the E-values: prec_mass x ppo_ (without
  // the minimum tolerance), or the Dalton tolerance when one is set.
  double compPrecErrorTole(double prec_mass) const;

  // consider zero ptm relaxed error
  double compRelaxErrorTole(double m1, double m2) const;

  double getPpo() const { return ppo_; }

  int getIntPpm() const;

  bool isUseMinTolerance() const { return use_min_tolerance_; }

  double getMinTolerance() const { return min_tolerance_; }

  void setPpo(double ppo) { ppo_ = ppo; }

  void setUseMinTolerance(bool use_min_tolerance) {
    use_min_tolerance_ = use_min_tolerance;
  }

  void setMinTolerance(double min_tolerance) { min_tolerance_ = min_tolerance; }

  // A tolerance in Dalton in place of mass x ppo_ (toppic's hidden
  // --prec-dalton-error, on the precursor tolerance only).
  void setDaltonTolerance(double tolerance) {
    use_dalton_tolerance_ = true;
    dalton_tolerance_ = tolerance;
  }

  bool isUseDaltonTolerance() const { return use_dalton_tolerance_; }

  double getDaltonTolerance() const { return dalton_tolerance_; }

 private:
  double ppo_;
  bool use_dalton_tolerance_ = false;
  double dalton_tolerance_ = 0;
  /* whether or not use minimum tolerance */
  bool use_min_tolerance_ = true;
  double min_tolerance_ = 0.01;
};

using PeakTolerancePtr = std::shared_ptr<PeakTolerance>;
using PeakTolerancePtrVec = std::vector<PeakTolerancePtr>;

}  // namespace toppic

#endif
