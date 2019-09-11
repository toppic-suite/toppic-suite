//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#ifndef TOPPIC_FILTER_MASS_MATCH_MASS_MATCH_FACTORY_HPP_
#define TOPPIC_FILTER_MASS_MATCH_MASS_MATCH_FACTORY_HPP_

#include "seq/proteoform.hpp"
#include "filter/massmatch/mass_match.hpp"

namespace toppic {

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

} /* namespace toppic */

#endif 
