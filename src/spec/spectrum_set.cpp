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
            = PrmMsFactory::geneMsTwoPtrVec(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
        srm_ms_two_ptr_vec_ 
            = PrmMsFactory::geneSuffixMsTwoPtrVec(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
        prm_ms_six_ptr_vec_ 
            = PrmMsFactory::geneMsSixPtrVec(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
      }
    }

PrmMsPtrVec SpectrumSet::getMsTwoPtrVec(SpParaPtr sp_para_ptr) {
  return PrmMsFactory::geneMsTwoPtrVec(deconv_ms_ptr_vec_,
                                       sp_para_ptr,
                                       prec_mono_mass_);
}

ExtendMsPtrVec SpectrumSet::getMsThreePtrVec(SpParaPtr sp_para_ptr) {
  return ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec_,
                                            sp_para_ptr,
                                            prec_mono_mass_ - std::accumulate(sp_para_ptr->mod_mass_.begin(), sp_para_ptr->mod_mass_.end(), 0.0));
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
  LOG_DEBUG(std::endl << "peak_num " << peak_num << std::endl);
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
