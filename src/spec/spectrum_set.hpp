/*
 * spectrum_set.hpp
 *
 *  Created on: Dec 9, 2013
 *      Author: xunlikun
 */

#ifndef PROT_SPECTRUM_SET_HPP_
#define PROT_SPECTRUM_SET_HPP_

#include <memory>
#include <vector>

#include "spec/deconv_ms.hpp"
#include "spec/sp_para.hpp"
#include "spec/extend_peak.hpp"
#include "spec/prm_peak.hpp"

namespace prot {

class SpectrumSet {
 public:
  SpectrumSet(DeconvMsPtr sp, double delta,SpParaPtr sp_para);
  DeconvMsPtr getDeconvMs(){return deconv_sp_;}
  ExtendMsPtr getSpThree(){return extend_ms_three_;}
  PrmMsPtr getSpTwo(){return prm_ms_two_;}
  PrmMsPtr getSpSix(){return prm_ms_six_;}
  PrmMsPtr getSpShiftSix(double shift){
    return getShiftSpSix(deconv_sp_,delta_,-shift,sp_para_ptr_);
  }
  double getDelta(){return delta_;}
 private:
  DeconvMsPtr deconv_sp_;
  SpParaPtr sp_para_ptr_;
  double delta_;
  PrmMsPtr prm_ms_two_;
  ExtendMsPtr extend_ms_three_;
  PrmMsPtr prm_ms_six_;
};

typedef std::shared_ptr<SpectrumSet> SpectrumSetPtr;

SpectrumSetPtr getSpectrumSet(DeconvMsPtr spectrum,double delta,
                              SpParaPtr sp_para);

} /* namespace prot */

#endif /* SPECTRUM_SET_HPP_ */
