#ifndef PROT_PRSM_PEAK_ION_PAIR_FACTORY_HPP_
#define PROT_PRSM_PEAK_ION_PAIR_FACTORY_HPP_

#include "base/proteoform.hpp"
#include "spec/extend_ms.hpp"
#include "spec/theo_peak.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

class PeakIonPairFactory {
 public:

  static PeakIonPairPtrVec findPairs(ExtendMsPtr ms_three_ptr, TheoPeakPtrVec &theo_peak_ptrs, 
                                  int bgn, int end, double add_tolerance);

  static PeakIonPairPtrVec genePeakIonPairs(const ProteoformPtr &proteoform_ptr, 
                                            const ExtendMsPtr &ms_three_ptr, double min_mass);

  static PeakIonPairPtrVec genePeakIonPairs(const ProteoformPtr &proteoform_ptr, 
                                            const ExtendMsPtrVec &ms_ptr_vec, double min_mass);
};

}
#endif
