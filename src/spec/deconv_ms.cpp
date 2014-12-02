
#include <memory>

#include "spec/peak.hpp"
#include "spec/ms.hpp"
#include "spec/deconv_ms.hpp"

namespace prot {


DeconvMsPtr getRefineMs(DeconvMsPtr deconv_ms_ptr, double calibration, 
                        double new_prec_mass) {
  MsHeaderPtr header_ptr = getHeaderPtr(deconv_ms_ptr, new_prec_mass);
  std::vector<DeconvPeakPtr> peak_ptr_list;
  for (size_t i = 0; i < deconv_ms_ptr->size(); i++) {
    DeconvPeakPtr ori_peak_ptr = deconv_ms_ptr->getPeakPtr(i);
    // * is a dereference operator
    DeconvPeakPtr new_peak_ptr(new DeconvPeak(*ori_peak_ptr.get()));
    new_peak_ptr->setPosition(new_peak_ptr->getPosition() * (1 + calibration));
    peak_ptr_list.push_back(new_peak_ptr);
  }
  DeconvMsPtr ms_ptr(new Ms<DeconvPeakPtr>(header_ptr, peak_ptr_list));
  return ms_ptr;
}

MsHeaderPtr getHeaderPtr(DeconvMsPtr deconv_ms_ptr, double new_prec_mass) {
  MsHeaderPtr header_ptr = deconv_ms_ptr->getHeaderPtr();
  MsHeaderPtr new_header_ptr(new MsHeader(*header_ptr.get()));
  double mono_mz = compMonoMz(new_prec_mass, header_ptr->getPrecCharge());
  new_header_ptr->setPrecMonoMz(mono_mz);
  return new_header_ptr;
}
    
/*
MsHeaderPtr getPrecHeaderPtr(DeconvMsPtr deconv_ms_ptr, double delta) {
  MsHeaderPtr header_ptr = deconv_ms_ptr->getHeaderPtr();
  MsHeaderPtr new_header_ptr(new MsHeader(*header_ptr.get()));
  double mono_mz = header_ptr->getPrecMonoMz() + delta
        / header_ptr->getPrecCharge();
  new_header_ptr->setPrecMonoMz(mono_mz);
  return new_header_ptr;
}
*/

}
