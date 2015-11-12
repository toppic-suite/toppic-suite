#ifndef PROT_SPEC_EXTEND_MS_FACTORY_HPP_
#define PROT_SPEC_EXTEND_MS_FACTORY_HPP_

#include "spec/deconv_ms.hpp"
#include "spec/extend_ms.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class ExtendMsFactory {

  static ExtendMsPtr geneMsThreePtr(DeconvMsPtr deconv_ms_ptr, SpParaPtr sp_para_ptr,
                                    double new_prec_mass);

  static ExtendMsPtrVec geneMsThreePtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                          SpParaPtr sp_para_ptr, double new_prec_mass);
};

} /* namespace prot */

#endif 
