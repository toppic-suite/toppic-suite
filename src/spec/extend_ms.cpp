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


#include <cmath>

#include "base/mass_constant.hpp"
#include "spec/extend_ms.hpp"

namespace prot {


std::vector<double> ExtendMs::getExtendMassVec (ExtendMsPtr extend_ms_ptr) {
  std::vector<double> masses;
  ExtendPeakPtrVec peak_ptr_list = extend_ms_ptr->getPeakPtrVec();
  for (size_t i = 0; i < peak_ptr_list.size(); i++) {
    masses.push_back(peak_ptr_list[i]->getPosition());
  }
  return masses;
}


inline bool massErrorUp(const std::pair<int, int> &a, const std::pair<int,int> b) {
    return a.first < b.first;
}

std::vector<std::pair<int, int>> ExtendMs::getExtendIntMassErrorList(
    const ExtendMsPtrVec &ext_ms_ptr_vec, bool pref, double scale){
  std::vector<std::pair<int,int>> mass_errors;
  for (size_t i = 0; i < ext_ms_ptr_vec.size(); i++) {
    ExtendMsPtr ext_ms_ptr = ext_ms_ptr_vec[i];
    double shift = 0;
    if (pref) {
      shift = ext_ms_ptr->getMsHeaderPtr()->getActivationPtr()->getNShift();
    }
    else {
      shift = ext_ms_ptr->getMsHeaderPtr()->getActivationPtr()->getCShift() 
          + MassConstant::getWaterMass();
    }
    std::pair<int,int> last_mass_error(-1, 0);
    for(size_t j=0; j<ext_ms_ptr->size(); j++){
      int m = (int) std::round((ext_ms_ptr->getPeakPtr(j)->getPosition() - shift) *scale);
      int e = (int) std::ceil(ext_ms_ptr->getPeakPtr(j)->getOrigTolerance()*scale);
      std::pair<int,int> cur_m_e (m, e);
      if(cur_m_e.first != last_mass_error.first){
        mass_errors.push_back(cur_m_e);
        last_mass_error = cur_m_e;
      }
      else if(cur_m_e.second > last_mass_error.second){
        mass_errors.pop_back();
        mass_errors.push_back(cur_m_e);
        last_mass_error = cur_m_e;
      }
    }
  }

  std::sort(mass_errors.begin(), mass_errors.end(),massErrorUp);
  return mass_errors;
}

} /* namespace prot */
