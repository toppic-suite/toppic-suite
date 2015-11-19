#ifndef PROT_PRSM_PEAK_ION_PAIR_FACTORY_HPP_
#define PROT_PRSM_PEAK_ION_PAIR_FACTORY_HPP_

#include "base/proteoform.hpp"
#include "spec/extend_ms.hpp"
#include "spec/theo_peak.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

class PeakIonPairFactory {
 public:
  PeakIonPairPtrVec getPeakIonPairs(const ProteoformPtr &proteoform_ptr, 
                                    const ExtendMsPtr &ms_three_ptr, double min_mass);

  PeakIonPairPtrVec getPeakIonPairs(const ProteoformPtr &proteoform_ptr, 
                                    const ExtendMsPtrVec &ms_ptr_vec, double min_mass);
};

}
#endif
