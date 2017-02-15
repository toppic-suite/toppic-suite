// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
