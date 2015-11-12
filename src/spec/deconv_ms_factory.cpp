#include <memory>

#include "spec/peak.hpp"
#include "spec/ms.hpp"
#include "spec/deconv_ms_factory.hpp"

namespace prot {


DeconvMsPtrVec DeconvMsFactory::getRefineMsPtrVec(const DeconvMsPtrVec &deconv_ms_ptr_vec, 
                                                  double new_prec_mass) {
  DeconvMsPtrVec result_ptrs;
  for (size_t m = 0; m < deconv_ms_ptr_vec.size(); m++) {
    DeconvMsPtr deconv_ms_ptr = deconv_ms_ptr_vec[m];
    MsHeaderPtr ori_header_ptr = deconv_ms_ptr->getMsHeaderPtr();
    MsHeaderPtr header_ptr = MsHeader::geneMsHeaderPtr(ori_header_ptr, new_prec_mass);
    std::vector<DeconvPeakPtr> peak_ptr_list;
    for (size_t p = 0; p < deconv_ms_ptr->size(); p++) {
      DeconvPeakPtr ori_peak_ptr = deconv_ms_ptr->getPeakPtr(p);
      // * is a dereference operator
      DeconvPeakPtr new_peak_ptr(new DeconvPeak(*ori_peak_ptr.get()));
      //new_peak_ptr->setPosition(new_peak_ptr->getPosition() * (1 + calibration));
      peak_ptr_list.push_back(new_peak_ptr);
    }
    DeconvMsPtr ms_ptr(new Ms<DeconvPeakPtr>(header_ptr, peak_ptr_list));
    result_ptrs.push_back(ms_ptr);
  }
  return result_ptrs;
}

}
