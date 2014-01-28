#ifndef PROT_DECONV_MS_HPP_
#define PROT_DECONV_MS_HPP_

#include <memory>

#include "spec/ms.hpp"
#include "spec/deconv_peak.hpp"

namespace prot {

typedef std::shared_ptr<Ms<DeconvPeakPtr>> DeconvMsPtr;

DeconvMsPtr  getRefineMs(DeconvMsPtr deconv_ms_ptr,
                         double calibration, double newPrecMass);

MsHeaderPtr getHeaderPtr(DeconvMsPtr deconv_ms_ptr,
                         double new_prec_mass);
    
MsHeaderPtr getDeltaHeaderPtr(DeconvMsPtr deconv_ms_ptr, double delta);

}

#endif
