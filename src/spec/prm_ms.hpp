#ifndef PROT_SPEC_PRM_MS_HPP_
#define PROT_SPEC_PRM_MS_HPP_

#include <memory>
#include <vector>

#include "spec/prm_peak.hpp"
#include "spec/ms.hpp"
#include "spec/peak_tolerance.hpp"

namespace prot {

typedef std::shared_ptr<Ms<PrmPeakPtr>> PrmMsPtr;
typedef std::vector<PrmMsPtr> PrmMsPtrVec;

class PrmMs {

  static std::vector<std::pair<int, int>> getIntMassErrorList(const PrmMsPtrVec &prm_ms_ptr_vec, 
                                                              PeakTolerancePtr tole_ptr,
                                                              double scale, bool n_strict, bool c_strict);

  static PrmPeakPtrVec getPrmPeakPtrs(const PrmMsPtrVec &prm_ms_ptr_vec, PeakTolerancePtr tole_ptr);
};

} /* namespace prot */

#endif 
