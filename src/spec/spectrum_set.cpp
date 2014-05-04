/*
 * spectrum_set.cpp
 *
 *  Created on: Dec 9, 2013
 *      Author: xunlikun
 */

#include "base/base_data.hpp"
#include "spec/extend_peak.hpp"
#include "spec/spectrum_set.hpp"

namespace prot {

SpectrumSet::SpectrumSet(const DeconvMsPtr &sp, double delta,
                         const SpParaPtr &sp_para_ptr){
  deconv_sp_ = sp;
  sp_para_ptr_ = sp_para_ptr;
  delta_=delta;
  prm_ms_two_ = getMsTwo(sp, delta, sp_para_ptr);
  extend_ms_three_ = getMsThree(sp,delta,sp_para_ptr);
  prm_ms_six_ = getMsSix(sp,delta,sp_para_ptr);
}

SpectrumSetPtr getSpectrumSet(const DeconvMsPtr &spectrum, double delta,
                              const SpParaPtr &sp_para_ptr){
  if((int)spectrum->size() < sp_para_ptr->getMinPeakNum() 
     || spectrum->getHeaderPtr()->getPrecMonoMass() < sp_para_ptr->getMinMass()){
    return SpectrumSetPtr(nullptr);
  }
  if(spectrum->getHeaderPtr()->getActivationPtr() == nullptr){
    if(sp_para_ptr->getActivationPtr()!=nullptr){
      spectrum->getHeaderPtr()->setActivationPtr(sp_para_ptr->getActivationPtr());
    }
    else{
      //logger
      return SpectrumSetPtr(nullptr);
    }
  }

  return SpectrumSetPtr(new SpectrumSet(spectrum,delta,sp_para_ptr));
}

} /* namespace prot */
