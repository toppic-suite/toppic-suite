#include "base/base_data.hpp"
#include "spec/extend_peak.hpp"
#include "spec/spectrum_set.hpp"

namespace prot {

SpectrumSet::SpectrumSet(DeconvMsPtr deconv_ms_ptr, double delta, 
                         SpParaPtr sp_para_ptr) {
  deconv_ms_ptr_ = deconv_ms_ptr;
  sp_para_ptr_ = sp_para_ptr;
  // add error tolerance for precursor mass 
  deconv_ms_ptr_->getHeaderPtr()->setErrorToleranceByPpo(
      sp_para_ptr_->getPeakTolerancePtr()->getPpo());
  delta_ = delta;
  prm_ms_two_ptr_ = createMsTwoPtr(deconv_ms_ptr, delta, sp_para_ptr);
  extend_ms_three_ptr_ = createMsThreePtr(deconv_ms_ptr,delta,sp_para_ptr);
  prm_ms_six_ptr_ = createMsSixPtr(deconv_ms_ptr,delta,sp_para_ptr);
}

SpectrumSetPtr getSpectrumSet(DeconvMsPtr deconv_ms_ptr, double delta,
                              SpParaPtr sp_para_ptr){
  if((int)deconv_ms_ptr->size() < sp_para_ptr->getMinPeakNum() 
     || deconv_ms_ptr->getHeaderPtr()->getPrecMonoMass() < sp_para_ptr->getMinMass()){
    return SpectrumSetPtr(nullptr);
  }
  if(deconv_ms_ptr->getHeaderPtr()->getActivationPtr() == nullptr){
    if(sp_para_ptr->getActivationPtr()!=nullptr){
      deconv_ms_ptr->getHeaderPtr()->setActivationPtr(sp_para_ptr->getActivationPtr());
    }
    else{
      //logger
      return SpectrumSetPtr(nullptr);
    }
  }

  return SpectrumSetPtr(new SpectrumSet(deconv_ms_ptr,delta,sp_para_ptr));
}

} /* namespace prot */
