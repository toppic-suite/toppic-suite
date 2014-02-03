/*
 * spectrum_set.cpp
 *
 *  Created on: Dec 9, 2013
 *      Author: xunlikun
 */

#include "spec/spectrum_set.hpp"
#include "base/base_data.hpp"

namespace prot {

SpectrumSet::SpectrumSet(DeconvMsPtr sp,double delta,SpParaPtr sp_para){
  deconv_sp_ = sp;
  sp_para_ptr_ = sp_para;
  delta_=delta;
  prm_ms_two_ = getMsTwo(sp,delta,sp_para);
  extend_ms_three_ = getMsThree(sp,delta,sp_para);
  prm_ms_six_ = getMsSix(sp,delta,sp_para);
}

SpectrumSetPtr getSpectrumSet(DeconvMsPtr spectrum,double delta,
                              SpParaPtr sp_para){
  if((int)spectrum->size() < sp_para->getMinPeakNum() 
     || spectrum->getHeaderPtr()->getPrecMonoMass() < sp_para->getMinMass()){
    //logger
    return nullptr;
  }
  if(spectrum->getHeaderPtr()->getActivationPtr() == nullptr){
    if(sp_para->getActivation()!=nullptr){
      spectrum->getHeaderPtr()->setActivationPtr(sp_para->getActivation());
    }
    else{
      //logger
      return nullptr;
    }
  }

  return SpectrumSetPtr(new SpectrumSet(spectrum,delta,sp_para));
}

} /* namespace prot */
