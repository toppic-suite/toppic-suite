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


#include <algorithm>
#include <cmath>
#include <iostream>

#include "base/logger.hpp"
#include "spec/prm_peak_factory.hpp"
#include "spec/prm_ms.hpp"

namespace prot {

inline bool massErrorUp(const std::pair<int, int> &a, const std::pair<int,int> b) {
  return a.first < b.first;
}

inline std::pair<int, int> getMassError(PrmPeakPtr peak_ptr, double scale, 
                                        bool n_strict, bool c_strict) {
  int m = (int)std::round(peak_ptr->getPosition()*scale);
  int e = 0;
  if(n_strict && c_strict){
    e = (int) std::ceil(peak_ptr->getStrictTolerance()*scale);
  }
  else if(n_strict && !c_strict){
    e = (int) std::ceil(peak_ptr->getNStrictCRelaxTolerance()*scale);
  }
  else if(!n_strict && c_strict){
    e = (int) std::ceil(peak_ptr->getNRelaxCStrictTolerance()*scale);
  }
  std::pair<int, int> mass_error(m, e);
  return mass_error;
}

std::vector<std::pair<int, int>> PrmMs::getIntMassErrorList(
    const PrmMsPtrVec &prm_ms_ptr_vec, PeakTolerancePtr tole_ptr,
    double scale, bool n_strict, bool c_strict){
  std::vector<std::pair<int,int>> mass_errors;
  for (size_t i = 0; i < prm_ms_ptr_vec.size(); i++) {
    PrmMsPtr prm_ms_ptr = prm_ms_ptr_vec[i];
    std::pair<int,int> last_mass_error(-1, 0);
    for(size_t j=0; j<prm_ms_ptr->size(); j++){
      std::pair<int, int> cur_m_e = getMassError(prm_ms_ptr->getPeakPtr(j), scale, n_strict, c_strict);
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
    //add zero mass for each spectrum to increase the score for zero mass
    double prec_mass = prm_ms_ptr_vec[i]->getMsHeaderPtr()->getPrecMonoMass();
    PrmPeakPtr zero_prm_ptr = PrmPeakFactory::getZeroPeakPtr(i, prec_mass, tole_ptr, 1);
    mass_errors.push_back(getMassError(zero_prm_ptr, scale, n_strict, c_strict));
    //add prec mass for each spectrum 
    PrmPeakPtr prec_prm_ptr = PrmPeakFactory::getPrecPeakPtr(i, prec_mass, tole_ptr, 1);
    mass_errors.push_back(getMassError(prec_prm_ptr, scale, true, true));
  }

  std::sort(mass_errors.begin(), mass_errors.end(),massErrorUp);
  return mass_errors;
}

PrmPeakPtrVec PrmMs::getPrmPeakPtrs(const PrmMsPtrVec &prm_ms_ptr_vec, 
                                    PeakTolerancePtr tole_ptr) {
  PrmPeakPtrVec peak_list;
  for (size_t i = 0; i < prm_ms_ptr_vec.size(); i++) {
    for(size_t j= 0;j<prm_ms_ptr_vec[i]->size() ;j++){
      peak_list.push_back(prm_ms_ptr_vec[i]->getPeakPtr(j));
    }
  }
  // add zero 
  double prec_mass = prm_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  // use spec_id = 0 and score = group_spec_num (size of prm_ms_ptr_vec)
  PrmPeakPtr zero_prm_ptr = PrmPeakFactory::getZeroPeakPtr(0, prec_mass, tole_ptr, prm_ms_ptr_vec.size());
  peak_list.push_back(zero_prm_ptr);
  // add prec mass  
  PrmPeakPtr prec_prm_ptr = PrmPeakFactory::getPrecPeakPtr(0, prec_mass, tole_ptr, prm_ms_ptr_vec.size());
  peak_list.push_back(prec_prm_ptr);
  std::sort(peak_list.begin(), peak_list.end(), PrmPeak::cmpPosInc);
  for (size_t i = 0; i < peak_list.size(); i++) {
    peak_list[i]->setPeakId(i);
  }
  return peak_list;
}

} /* namespace prot */
