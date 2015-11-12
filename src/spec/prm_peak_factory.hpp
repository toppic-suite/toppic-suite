#ifndef PROT_SPEC_PRM_PEAK_FACTORY_HPP_
#define PROT_SPEC_PRM_PEAK_FACTORY_HPP_

#include "spec/peak_tolerance.hpp"
#include "spec/prm_peak.hpp"

namespace prot {

class PrmPeakFactory {
 public:
  static PrmPeakPtr getZeroPeakPtr(int spec_id, double prec_mono_mass, 
                                   PeakTolerancePtr tole_ptr, double score);

  static PrmPeakPtr getPrecPeakPtr(int spec_id, double prec_mono_mass, 
                                   PeakTolerancePtr tole_ptr, double score);
};

} /* namespace prot */

#endif 
