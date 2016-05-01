#ifndef ZERO_PTM_FILTER_MASS_MATCH_FACTORY_HPP_
#define ZERO_PTM_FILTER_MASS_MATCH_FACTORY_HPP_

#include <cmath>

#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "zeroptmfilter/mass_match.hpp"

namespace prot {

class MassMatchFactory {
 public:

  static MassMatchPtr getPrmDiagMassMatchPtr(const ProteoformPtrVec &proteo_ptrs,
                                             double max_proteoform_mass, double scale);

  static MassMatchPtr getSrmDiagMassMatchPtr(const ProteoformPtrVec &proteo_ptrs,
                                             std::vector<std::vector<double>> &n_ace_shift_2d,
                                             double max_proteoform_mass, double scale);

  static MassMatchPtr getPrmTermMassMatchPtr(const ProteoformPtrVec &proteo_ptrs, 
                                             std::vector<std::vector<double>> &real_shift_2d,
                                             double max_proteoform_mass, double scale);


  static MassMatchPtr getSrmTermMassMatchPtr(const ProteoformPtrVec &proteo_ptrs, 
                                             std::vector<std::vector<double>> &real_shift_2d,
                                             std::vector<std::vector<double>> &n_ace_shift_2d,
                                             double max_proteoform_mass, double scale);
};

} /* namespace prot */

#endif 
