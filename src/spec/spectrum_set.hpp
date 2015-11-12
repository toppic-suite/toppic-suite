#ifndef PROT_SPEC_SPECTRUM_SET_HPP_
#define PROT_SPEC_SPECTRUM_SET_HPP_

#include <memory>
#include <vector>

#include "spec/deconv_ms.hpp"
#include "spec/extend_ms.hpp"
#include "spec/prm_ms.hpp"
#include "spec/prm_ms_factory.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class SpectrumSet;
typedef std::shared_ptr<SpectrumSet> SpectrumSetPtr;

class SpectrumSet {
 public:
  SpectrumSet(DeconvMsPtrVec deconv_ms_ptr_vec, SpParaPtr sp_para_ptr, 
              double prec_mono_mass);

  double getPrecMonoMass() {return prec_mono_mass_;}

  bool isValid() {return valid_;}

  int getSpecId() {return deconv_ms_ptr_vec_[0]->getMsHeaderPtr()->getId();}

  ExtendMsPtrVec getMsThreePtrVec() {return extend_ms_three_ptr_vec_;}

  DeconvMsPtrVec getDeconvMsPtrVec(){return deconv_ms_ptr_vec_;}

  PrmMsPtrVec getMsTwoPtrVec() {return prm_ms_two_ptr_vec_;}

  PrmMsPtrVec getMsSixPtrVec(){return prm_ms_six_ptr_vec_;}

  PrmMsPtrVec getMsShiftSixPtrVec(double shift){
    return PrmMsFactory::geneShiftMsSixPtrVec(
        deconv_ms_ptr_vec_, sp_para_ptr_, prec_mono_mass_, -shift);
  }

 private:
  DeconvMsPtrVec deconv_ms_ptr_vec_;
  SpParaPtr sp_para_ptr_;
  double prec_mono_mass_;
  bool valid_ = true;
  ExtendMsPtrVec extend_ms_three_ptr_vec_;
  PrmMsPtrVec prm_ms_two_ptr_vec_;
  PrmMsPtrVec prm_ms_six_ptr_vec_;

  bool checkValid(SpParaPtr sp_para_ptr);
};

} /* namespace prot */

#endif /* SPECTRUM_SET_HPP_ */
