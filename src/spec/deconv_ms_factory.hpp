#ifndef PROT_SPEC_DECONV_MS_FACTORY_HPP_
#define PROT_SPEC_DECONV_MS_FACTORY_HPP_

#include "spec/ms.hpp"
#include "spec/deconv_peak.hpp"

namespace prot {

typedef std::shared_ptr<Ms<DeconvPeakPtr>> DeconvMsPtr;

typedef std::vector<DeconvMsPtr> DeconvMsPtrVec;

class DeconvMsFactory {
 public:
  static DeconvMsPtrVec  getRefineMsPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                           double new_prec_mass);
};

}

#endif
