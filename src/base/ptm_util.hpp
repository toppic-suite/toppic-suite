#ifndef PROT_PTM_UTIL_H_
#define PROT_PTM_UTIL_H_

#include <vector>

#include "ptm.hpp"

namespace prot {


PtmPtrVec getPtmPtrVecInstance(const char* file_name);
/**
 * Returns a PTM based on the abbreviation name. Returns null if the
 * abbreviation name does not exist.
 */
PtmPtr getPtmPtrByAbbrName(PtmPtrVec &ptm_ptr_vec, 
                           const std::string &abbr_name);

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool containAbbrsName(PtmPtrVec &ptm_ptr_vec, const std::string &abbr_name);

PtmPtr findEmptyPtmPtr(PtmPtrVec &ptm_ptr_vec);


}
#endif
