#ifndef ZERO_PTM_FILTER_MASS_MATCH_FACTORY_HPP_
#define ZERO_PTM_FILTER_MASS_MATCH_FACTORY_HPP_

#include <cmath>

#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "zeroptmfilter/mass_match.hpp"

namespace prot {

class MassMatchFactory {
 public:
  static MassMatchPtr getMassMatchPtr(const ProteoformPtrVec &proteo_ptrs,  
                                      double max_proteoform_mass, double scale, bool rev);  

  static MassMatchPtr getMassMatchPtr(const ProteoformPtrVec &proteo_ptrs, 
                                      std::vector<std::vector<double>> &shift_2d,
                                      double max_proteoform_mass, double scale, bool rev);
};

} /* namespace prot */

#endif 
