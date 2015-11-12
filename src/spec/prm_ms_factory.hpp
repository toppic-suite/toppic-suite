#ifndef PROT_SPEC_PRM_MS_FACTORY_HPP_
#define PROT_SPEC_PRM_MS_FACTORY_HPP_

#include "spec/sp_para.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/prm_ms.hpp"

namespace prot {

class PrmMsFactory {
 public:
  static PrmMsPtrVec geneMsTwoPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                     SpParaPtr sp_para_ptr,
                                     double prec_mono_mass);

  static PrmMsPtrVec geneMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                     SpParaPtr sp_para_ptr,
                                     double prec_mono_mass);

  static PrmMsPtrVec geneShiftMsSixPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                          SpParaPtr sp_para_ptr, 
                                          double prec_mono_mass, double shift);
};

} /* namespace prot */

#endif /* PRM_PEAK_HPP_ */
