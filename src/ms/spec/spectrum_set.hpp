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

#ifndef TOPPIC_MS_SPEC_SPECTRUM_SET_HPP_
#define TOPPIC_MS_SPEC_SPECTRUM_SET_HPP_

#include <memory>
#include <vector>

#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/extend_ms.hpp"
#include "ms/spec/prm_ms.hpp"

namespace toppic {

class SpectrumSet;
using SpectrumSetPtr = std::shared_ptr<SpectrumSet>;
using SpectrumSetPtrVec = std::vector<SpectrumSetPtr>;

class SpectrumSet {
 public:
  SpectrumSet(const DeconvMsPtrVec& deconv_ms_ptr_vec, double prec_mono_mass,
              double n_term_label_mass, bool valid,
              const ExtendMsPtrVec& extend_ms_three_ptr_vec,
              const PrmMsPtrVec& prm_ms_two_ptr_vec,
              const PrmMsPtrVec& srm_ms_two_ptr_vec,
              const PrmMsPtrVec& prm_ms_six_ptr_vec);

  double getPrecMonoMass() const { return prec_mono_mass_; }

  double getNTermLabelMass() const { return n_term_label_mass_; }

  bool isValid() const { return valid_; }

  int getSpectrumId() const {
    return deconv_ms_ptr_vec_[0]->getMsHeaderPtr()->getSpecId();
  }

  const ExtendMsPtrVec& getMsThreePtrVec() const {
    return extend_ms_three_ptr_vec_;
  }

  const DeconvMsPtrVec& getDeconvMsPtrVec() const { return deconv_ms_ptr_vec_; }

  const PrmMsPtrVec& getMsTwoPtrVec() const { return prm_ms_two_ptr_vec_; }

  const PrmMsPtrVec& getSuffixMsTwoPtrVec() const {
    return srm_ms_two_ptr_vec_;
  }

  const PrmMsPtrVec& getMsSixPtrVec() const { return prm_ms_six_ptr_vec_; }

 private:
  DeconvMsPtrVec deconv_ms_ptr_vec_;

  double prec_mono_mass_;

  double n_term_label_mass_ = 0;

  bool valid_ = true;

  ExtendMsPtrVec extend_ms_three_ptr_vec_;

  PrmMsPtrVec prm_ms_two_ptr_vec_;

  PrmMsPtrVec srm_ms_two_ptr_vec_;

  PrmMsPtrVec prm_ms_six_ptr_vec_;
};

}  // namespace toppic

#endif
