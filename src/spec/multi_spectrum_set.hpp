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
  MultiSpectrumSet(DeconvMsPtr deconv_ms_ptr, double delta, SpParaPtr sp_para_ptr);

  DeconvMsPtrVec getDeconvMsPtrVec(){return deconv_sp_vec_;}

  ExtendMsPtr getSpThreePtr(){return extend_ms_three_ptr_;}

  PrmMsPtr getSpTwoPtr(){return prm_ms_two_ptr_;}

  PrmMsPtr getSpSixPtr(){return prm_ms_six_ptr_;}

  PrmMsPtr getSpShiftSixPtr(double shift){
    return getShiftMsSixPtr(deconv_ms_ptr_, delta_, -shift, sp_para_ptr_);
  }
  double getDelta(){return delta_;}

 private:
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  double adjusted_prec_mass_;
  SpParaPtr sp_para_ptr_;
  PrmMsPtr prm_ms_two_ptr_;
  ExtendMsPtr extend_ms_three_ptr_;
  PrmMsPtr prm_ms_six_ptr_;
};

typedef std::shared_ptr<SpectrumSet> SpectrumSetPtr;

SpectrumSetPtr getSpectrumSet(DeconvMsPtr deconv_ms_ptr, double delta,
                              SpParaPtr sp_para_ptr);

} /* namespace prot */

#endif /* SPECTRUM_SET_HPP_ */
