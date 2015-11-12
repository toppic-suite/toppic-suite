#include <functional>

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
      double ppo = sp_para_ptr_->getPeakTolerancePtr()->getPpo();
      ActivationPtr activation_ptr = sp_para_ptr->getActivationPtr();
      for (size_t i = 0; i < deconv_ms_ptr_vec_.size(); i++) {
        DeconvMsPtr deconv_ms_ptr = deconv_ms_ptr_vec[i];
        MsHeaderPtr ms_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
        ms_header_ptr->setErrorToleranceByPpo(ppo);
        if(ms_header_ptr->getActivationPtr() == nullptr && activation_ptr != nullptr){
          ms_header_ptr->setActivationPtr(activation_ptr);
        }
      }
      valid_ = checkValid(sp_para_ptr);
      if (valid_) {
        extend_ms_three_ptr_vec_ 
            = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec_, sp_para_ptr, prec_mono_mass);
        prm_ms_two_ptr_vec_ 
            = PrmMsFactory::geneMsTwoPtrVec(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
        prm_ms_six_ptr_vec_ 
            = PrmMsFactory::geneMsSixPtrVec(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
      }
    }

bool SpectrumSet::checkValid(SpParaPtr sp_para_ptr) {
  if (prec_mono_mass_ < sp_para_ptr->getMinMass()) {
    return false;
  }
  int peak_num = 0;
  for (size_t i = 0; i < deconv_ms_ptr_vec_.size(); i++) {
    peak_num += deconv_ms_ptr_vec_[i]->size();
  }
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

} /* namespace prot */
