/*
 * spectrum_set.hpp
 *
 *  Created on: Dec 9, 2013
 *      Author: xunlikun
 */

#ifndef PROT_MULTI_SPECTRUM_SET_HPP_
#define PROT_MULTI_SPECTRUM_SET_HPP_

#include <memory>
#include <vector>

#include "spec/deconv_ms.hpp"
#include "spec/sp_para.hpp"
#include "spec/extend_peak.hpp"
#include "spec/prm_peak.hpp"

namespace prot {

class MultiSpectrumSet {
 public:
  MultiSpectrumSet(const DeconvMsPtr &sp, double delta,
                   const SpParaPtr &sp_para);
  DeconvMsPtrVec getDeconvMsPtrVec(){return deconv_sp_vec_;}
  ExtendMsPtr getSpThree(){return extend_ms_three_;}
  PrmMsPtr getSpTwo(){return prm_ms_two_;}
  PrmMsPtr getSpSix(){return prm_ms_six_;}
  PrmMsPtr getSpShiftSix(double shift){
    return getShiftMsSix(deconv_sp_,delta_,-shift,sp_para_ptr_);
  }
  double getDelta(){return delta_;}

 private:
  DeconvMsPtrVec deconv_sp_vec_;
  double adjusted_prec_mass;
  SpParaPtr sp_para_ptr_;
  PrmMsPtr prm_ms_two_;
  ExtendMsPtr extend_ms_three_;
  PrmMsPtr prm_ms_six_;
};

typedef std::shared_ptr<SpectrumSet> SpectrumSetPtr;

SpectrumSetPtr getSpectrumSet(const DeconvMsPtr &spectrum,
                              double delta,
                              const SpParaPtr &sp_para);

} /* namespace prot */

#endif /* SPECTRUM_SET_HPP_ */
