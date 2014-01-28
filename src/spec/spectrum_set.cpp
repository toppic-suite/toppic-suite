/*
 * spectrum_set.cpp
 *
 *  Created on: Dec 9, 2013
 *      Author: xunlikun
 */

#include <spec/spectrum_set.hpp>
#include "base/base_data.hpp"

namespace prot {

SpectrumSet::SpectrumSet(DeconvMsPtr sp,double delta,SpParaPtr sp_para,
                         double shift){
    deconv_sp_ = sp;
    delta_=delta;
    prm_ms_two_ = prot::getMsTwo(sp,delta,sp_para);
    extend_ms_three_ = prot::getMsThree(sp,delta,sp_para);
    prm_ms_six_=prot::getSpSix(sp,delta,sp_para);
    prm_ms_shift_six_=prot::getShiftSpSix(sp,delta,-shift,sp_para);
}

SpectrumSetPtr getSpectrumSet(DeconvMsPtr spectrum,double delta,
                              SpParaPtr sp_para,double shift){
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

    return SpectrumSetPtr(new SpectrumSet(spectrum,delta,sp_para,shift));
}

} /* namespace prot */
