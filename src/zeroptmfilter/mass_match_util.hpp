#ifndef ZERO_PTM_FILTER_MASS_MATCH_UTIL_HPP_
#define ZERO_PTM_FILTER_MASS_MATCH_UTIL_HPP_

#include <cmath>

#include "zeroptmfilter/filter_protein.hpp"
#include "zeroptmfilter/mass_match.hpp"

namespace prot {

class MassMatchUtil {
 public:
  static FilterProteinPtrVec findTopProteins(std::vector<short> &scores, 
                                             std::vector<int> &proteo_row_begins,
                                             std::vector<int> &proteo_row_ends,
                                             int threshold, int num);

  static FilterProteinPtrVec findTopProteins(std::vector<short> &scores, 
                                             std::vector<short> &rev_scores, 
                                             MassMatchPtr index_ptr,
                                             MassMatchPtr rev_index_ptr,
                                             double threshold, int num,
                                             bool add_shifts, int shift_num);
};

} /* namespace prot */

#endif 
