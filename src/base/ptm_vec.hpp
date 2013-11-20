#ifndef PROTOMICS_PTM_VEC_H_
#define PROTOMICS_PTM_VEC_H_

#include <vector>

#include "ptm.hpp"

namespace proteomics {


class PtmVec {
 public:
  static const std::vector<PtmPtr> getInstance(std::vector<AcidPtr> &acid_ptrs, 
                                               const char* file_name);
  /**
   * Returns a PTM based on the abbreviation name. Returns null if the
   * abbreviation name does not exist.
   */
  static PtmPtr getPtmPtrByAbbrName(std::vector<PtmPtr> &ptm_ptrs, 
                                    const std::string &abbr_name);

  /**
   * Checks if the list contains an amino acid with the specific name.
   */
  static bool containAbbrsName(std::vector<PtmPtr> &ptm_ptrs, const std::string &abbr_name);

  static PtmPtr getEmptyPtmPtr(std::vector<PtmPtr> &ptm_ptrs);

  static std::vector<AcidPtr> getValidAcidPtrs(std::vector<PtmPtr> &ptm_ptrs);
};
}
#endif
