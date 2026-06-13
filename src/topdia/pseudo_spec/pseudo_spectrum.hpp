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

#ifndef TOPPIC_TOPDIA_PSEUDO_SPEC_PSEUDO_SPECTRUM_HPP_
#define TOPPIC_TOPDIA_PSEUDO_SPEC_PSEUDO_SPECTRUM_HPP_

#include <memory>
#include <utility>
#include <vector>

#include "topdia/pseudo_spec/mzrt_feature.hpp"
#include "topdia/pseudo_spec/pseudo_peak.hpp"

namespace toppic {

class PseudoSpectrum {
 public:
  PseudoSpectrum(MzrtFeaturePtr ms1Feature,
                 const std::vector<PseudoPeak>& fragmentFeatures)
      : ms1_feature_(std::move(ms1Feature)),
        fragment_features_(fragmentFeatures) {}

  void add_fragments(const std::vector<PseudoPeak>& fragmentFeatures);
  void add_fragment(const PseudoPeak& fragmentFeature);
  const MzrtFeaturePtr& getMs1Feature() const { return ms1_feature_; }
  std::vector<PseudoPeak>& getFragmentFeatures() { return fragment_features_; }
  int getAssignedFragmentCount() const { return fragment_features_.size(); }

 private:
  MzrtFeaturePtr ms1_feature_;
  std::vector<PseudoPeak> fragment_features_;
};
using PseudoSpectrumPtr = std::shared_ptr<PseudoSpectrum>;
using PseudoSpectrumPtrVec = std::vector<PseudoSpectrumPtr>;

}  // namespace toppic

#endif  // TOPPIC_TOPDIA_PSEUDO_SPEC_PSEUDO_SPECTRUM_HPP_
