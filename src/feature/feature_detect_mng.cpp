// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "base/mass_constant.hpp"
#include "feature/feature_detect_mng.hpp"

namespace prot {

FeatureDetectMng::FeatureDetectMng() {
  double ppo = 0.000015;
  bool use_min_tolerance = true;
  double min_tolerance = 0.01;
  peak_tolerance_ptr_ = PeakTolerancePtr(
      new PeakTolerance(ppo, use_min_tolerance, min_tolerance));

  // extend sp parameter 
  double IM = MassConstant::getIsotopeMass();
  // the set of offsets used to expand the monoisotopic mass list 
  std::vector<double> offsets {{0, -IM, IM, -2 * IM, 2 * IM}};
  ext_offsets_ = offsets;
}

std::vector<double> FeatureDetectMng::getExtMasses(double mass) {
  std::vector<double> result;
  if (mass < extend_min_mass_) {
    result.push_back(mass);
  }
  else {
    for (size_t i = 0; i < ext_offsets_.size(); i++) {
      double new_mass = mass + ext_offsets_[i];
      result.push_back(new_mass);
    }
  }
  return result;
}

} /* namespace */

