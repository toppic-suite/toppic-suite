//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#ifndef TOPPIC_MS_SPEC_SPECTRUM_SET_HPP_
#define TOPPIC_MS_SPEC_SPECTRUM_SET_HPP_

#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/extend_ms.hpp"
#include "ms/spec/prm_ms.hpp"

namespace toppic {

class SpectrumSet;
typedef std::shared_ptr<SpectrumSet> SpectrumSetPtr;
typedef std::vector<SpectrumSetPtr> SpectrumSetPtrVec;

class SpectrumSet {
 public:
  SpectrumSet(DeconvMsPtrVec deconv_ms_ptr_vec, 
              double prec_mono_mass,
              bool valid, 
              ExtendMsPtrVec extend_ms_three_ptr_vec,
              PrmMsPtrVec prm_ms_two_ptr_vec,
              PrmMsPtrVec srm_ms_two_ptr_vec,
              PrmMsPtrVec prm_ms_six_ptr_vec);

  double getPrecMonoMass() {return prec_mono_mass_;}

  bool isValid() {return valid_;}

  int getSpectrumId() {return deconv_ms_ptr_vec_[0]->getMsHeaderPtr()->getId();}

  ExtendMsPtrVec getMsThreePtrVec() {return extend_ms_three_ptr_vec_;}

  DeconvMsPtrVec getDeconvMsPtrVec() {return deconv_ms_ptr_vec_;}

  PrmMsPtrVec getMsTwoPtrVec() {return prm_ms_two_ptr_vec_;}

  PrmMsPtrVec getSuffixMsTwoPtrVec() {return srm_ms_two_ptr_vec_;}

  PrmMsPtrVec getMsSixPtrVec() {return prm_ms_six_ptr_vec_;}

  //PrmMsPtrVec getSuffixMsSixPtrVec() {return srm_ms_six_ptr_vec_;}

 private:
  DeconvMsPtrVec deconv_ms_ptr_vec_;

  double prec_mono_mass_;

  bool valid_ = true;

  ExtendMsPtrVec extend_ms_three_ptr_vec_;

  PrmMsPtrVec prm_ms_two_ptr_vec_;

  PrmMsPtrVec srm_ms_two_ptr_vec_;

  PrmMsPtrVec prm_ms_six_ptr_vec_;

  //suffix ms six is not used
  //PrmMsPtrVec srm_ms_six_ptr_vec_;

};

} /* namespace toppic */

#endif /* SPECTRUM_SET_HPP_ */
