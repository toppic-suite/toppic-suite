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
  SpectrumSet(DeconvMsPtrVec deconv_ms_ptr_vec, SpParaPtr sp_para_ptr, 
              double prec_mono_mass);

  double getPrecMonoMass() {return prec_mono_mass_;}

  bool isFiltered() {return filtered;}

  ExtendMsPtrVec getMsThreePtrVec() {return extend_ms_three_ptr_vec_;}

  DeconvMsPtrVec getDeconvMsPtrVec(){return deconv_ms_ptr_vec_;}

  /*
  PrmMsPtr getMsTwoPtr() {return prm_ms_two_ptr_;}

  PrmMsPtr getMsSixPtr(){return prm_ms_six_ptr_;}

  PrmMsPtr getMsShiftSixPtr(double shift){
    return createShiftMsSixPtr(deconv_ms_ptr_, delta_, -shift, sp_para_ptr_);
  }
  */


 private:
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  SpParaPtr sp_para_ptr_;
  double prec_mono_mass_;
  bool filtered = false;
  ExtendMsPtrVec extend_ms_three_ptr_vec_;

  bool checkFiltration(SpParaPtr sp_para_ptr);
  //PrmMsPtr prm_ms_two_ptr_;
  //PrmMsPtr prm_ms_six_ptr_;
};

typedef std::shared_ptr<SpectrumSet> SpectrumSetPtr;

} /* namespace prot */

#endif /* SPECTRUM_SET_HPP_ */
