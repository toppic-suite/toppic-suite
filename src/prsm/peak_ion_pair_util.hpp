#ifndef PROT_PRSM_PEAK_ION_PAIR_UTIL_HPP_
#define PROT_PRSM_PEAK_ION_PAIR_UTIL_HPP_

#include "spec/rm_break_type.hpp"
#include "spec/extend_peak.hpp"
#include "spec/theo_peak.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

class PeakIonPairUtil {
 public:
  static PeakIonPairPtrVec getMatchedPairs(const PeakIonPairPtrVec &pair_ptrs, 
                                           int spec_id, int peak_id);

  static int getPeakIonPairNum(PeakIonPairPtrVec pair_ptrs); 

  static double computePairConverage(const PeakIonPairPtrVec &pair_ptrs, int begin, 
                                     int end, RmBreakTypePtr type_ptr);

};

}
#endif
