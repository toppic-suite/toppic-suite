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
  SpectrumSet(DeconvMsPtr deconv_ms_ptr, double delta, SpParaPtr sp_para_ptr);

  DeconvMsPtr getDeconvMsPtr(){return deconv_ms_ptr_;}

  ExtendMsPtr getMsThreePtr() {return extend_ms_three_ptr_;}

  PrmMsPtr getMsTwoPtr() {return prm_ms_two_ptr_;}

  PrmMsPtr getMsSixPtr(){return prm_ms_six_ptr_;}

  PrmMsPtr getMsShiftSixPtr(double shift){
    return createShiftMsSixPtr(deconv_ms_ptr_, delta_, -shift, sp_para_ptr_);
  }

  double getDelta(){return delta_;}

 private:
  DeconvMsPtr deconv_ms_ptr_;
  SpParaPtr sp_para_ptr_;
  double delta_;
  PrmMsPtr prm_ms_two_ptr_;
  ExtendMsPtr extend_ms_three_ptr_;
  PrmMsPtr prm_ms_six_ptr_;
};

typedef std::shared_ptr<SpectrumSet> SpectrumSetPtr;

SpectrumSetPtr getSpectrumSet(DeconvMsPtr deconv_ms_ptr, double delta, 
                              SpParaPtr sp_para_ptr);

} /* namespace prot */

#endif /* SPECTRUM_SET_HPP_ */
