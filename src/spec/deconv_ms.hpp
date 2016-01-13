#ifndef PROT_SPEC_DECONV_MS_HPP_
#define PROT_SPEC_DECONV_MS_HPP_

#include "spec/ms.hpp"
#include "spec/deconv_peak.hpp"

namespace prot {

typedef std::shared_ptr<Ms<DeconvPeakPtr>> DeconvMsPtr;

typedef std::vector<DeconvMsPtr> DeconvMsPtrVec;

}

#endif
