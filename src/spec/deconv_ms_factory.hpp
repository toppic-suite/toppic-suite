#ifndef PROT_SPEC_DECONV_MS_FACTORY_HPP_
#define PROT_SPEC_DECONV_MS_FACTORY_HPP_

#include "spec/deconv_ms.hpp"

namespace prot {

class DeconvMsFactory {
 public:
  static DeconvMsPtrVec  getRefineMsPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                           double new_prec_mass);
};

}

#endif
