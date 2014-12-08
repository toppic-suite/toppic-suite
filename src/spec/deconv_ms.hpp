#ifndef PROT_DECONV_MS_HPP_
#define PROT_DECONV_MS_HPP_

#include "spec/ms.hpp"
#include "spec/deconv_peak.hpp"

namespace prot {

typedef std::shared_ptr<Ms<DeconvPeakPtr>> DeconvMsPtr;

typedef std::vector<DeconvMsPtr> DeconvMsPtrVec;

//DeconvMsPtr  getRefineMs(DeconvMsPtr deconv_ms_ptr, double new_prec_mass);

DeconvMsPtrVec  getRefineMsPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                  double new_prec_mass);

MsHeaderPtr getHeaderPtr(DeconvMsPtr deconv_ms_ptr, double new_prec_mass);
    
}

#endif
