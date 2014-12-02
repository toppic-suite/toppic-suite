#include "base/base_data.hpp"
#include "spec/extend_peak.hpp"
#include "spec/spectrum_set.hpp"

namespace prot {

SpectrumSet::SpectrumSet(DeconvMsPtrVec deconv_ms_ptr_vec, SpParaPtr sp_para_ptr,
                         double prec_mono_mass) {
  deconv_ms_ptr_vec_ = deconv_ms_ptr_vec;
  sp_para_ptr_ = sp_para_ptr;
  // add error tolerance for precursor mass 
  double ppo = sp_para_ptr_->getPeakTolerancePtr()->getPpo();
  ActivationPtr activation_ptr = sp_para_ptr->getActivationPtr();
  for (size_t i = 0; i < deconv_ms_ptr_vec_.size(); i++) {
    DeconvMsPtr deconv_ms_ptr = deconv_ms_ptr_vec[i];
    deconv_ms_ptr->getHeaderPtr()->setErrorToleranceByPpo(ppo);
    if(deconv_ms_ptr->getHeaderPtr()->getActivationPtr() == nullptr 
       && activation_ptr != nullptr){
      deconv_ms_ptr->getHeaderPtr()->setActivationPtr(activation_ptr);
    }
  }
  prec_mono_mass_ = prec_mono_mass;
  filtered = checkFiltration(sp_para_ptr);
  if (!filtered) {
    extend_ms_three_ptr_vec_ 
        = createMsThreePtrVec(deconv_ms_ptr_vec_, sp_para_ptr, prec_mono_mass);
  //prm_ms_two_ptr_ = createMsTwoPtr(deconv_ms_ptr, delta, sp_para_ptr);
  //prm_ms_six_ptr_ = createMsSixPtr(deconv_ms_ptr,delta,sp_para_ptr);
  }
}

bool SpectrumSet::checkFiltration(SpParaPtr sp_para_ptr) {
  if (prec_mono_mass_ < sp_para_ptr->getMinMass()) {
    return true;
  }
  int peak_num = 0;
  for (size_t i = 0; i < deconv_ms_ptr_vec_.size(); i++) {
    peak_num += deconv_ms_ptr_vec_[i]->size();
  }
  if(peak_num < sp_para_ptr->getMinPeakNum()){
    return true;
  }
  for (size_t i = 0; i < deconv_ms_ptr_vec_.size(); i++) {
    if(deconv_ms_ptr_vec_[i]->getHeaderPtr()->getActivationPtr() == nullptr){
      return true;
    }
  }
  return false;
}

} /* namespace prot */
