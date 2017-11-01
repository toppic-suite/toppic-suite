//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include <numeric>
#include <functional>

#include "base/logger.hpp"
#include "spec/extend_ms_factory.hpp"
#include "spec/prm_ms_factory.hpp"
#include "spec/spectrum_set.hpp"

namespace prot {

SpectrumSet::SpectrumSet(DeconvMsPtrVec deconv_ms_ptr_vec,
                         SpParaPtr sp_para_ptr,
                         double prec_mono_mass): 
    deconv_ms_ptr_vec_(deconv_ms_ptr_vec),
    sp_para_ptr_(sp_para_ptr),
    prec_mono_mass_(prec_mono_mass) {
      // add error tolerance for precursor mass 
      ActivationPtr activation_ptr = sp_para_ptr->getActivationPtr();
      for (size_t i = 0; i < deconv_ms_ptr_vec_.size(); i++) {
        DeconvMsPtr deconv_ms_ptr = deconv_ms_ptr_vec[i];
        MsHeaderPtr ms_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
        if(ms_header_ptr->getActivationPtr() == nullptr && activation_ptr != nullptr){
          ms_header_ptr->setActivationPtr(activation_ptr);
        }
      }
      valid_ = checkValid(sp_para_ptr);
      //LOG_DEBUG("valid " << valid_);
      if (valid_) {
        extend_ms_three_ptr_vec_ 
            = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec_, sp_para_ptr, prec_mono_mass);
        prm_ms_two_ptr_vec_ 
            = PrmMsFactory::geneMsTwoPtrVec(deconv_ms_ptr_vec_, sp_para_ptr, prec_mono_mass);
        srm_ms_two_ptr_vec_ 
            = PrmMsFactory::geneSuffixMsTwoPtrVec(deconv_ms_ptr_vec_, sp_para_ptr, prec_mono_mass);
        prm_ms_six_ptr_vec_ 
            = PrmMsFactory::geneMsSixPtrVec(deconv_ms_ptr_vec_, sp_para_ptr, prec_mono_mass);
      }
    }

PrmMsPtrVec SpectrumSet::getMsTwoPtrVec(SpParaPtr sp_para_ptr) {
  return PrmMsFactory::geneMsTwoPtrVec(deconv_ms_ptr_vec_,
                                       sp_para_ptr,
                                       prec_mono_mass_);
}

PrmMsPtrVec SpectrumSet::getSuffixMsTwoPtrVec(SpParaPtr sp_para_ptr) {
  return PrmMsFactory::geneSuffixMsTwoPtrVec(deconv_ms_ptr_vec_, sp_para_ptr, prec_mono_mass_);
}

bool SpectrumSet::checkValid(SpParaPtr sp_para_ptr) {
  if (prec_mono_mass_ < sp_para_ptr->getMinMass()) {
    return false;
  }
  int peak_num = 0;
  for (size_t i = 0; i < deconv_ms_ptr_vec_.size(); i++) {
    peak_num += deconv_ms_ptr_vec_[i]->size();
  }
  LOG_DEBUG("peak_num " << peak_num);
  if(peak_num < sp_para_ptr->getMinPeakNum()){
    return false;
  }
  for (size_t i = 0; i < deconv_ms_ptr_vec_.size(); i++) {
    if(deconv_ms_ptr_vec_[i]->getMsHeaderPtr()->getActivationPtr() == nullptr){
      return false;
    }
  }
  return true;
}

}  // namespace prot
