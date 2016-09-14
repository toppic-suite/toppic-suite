// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#include <iostream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/prot_mod.hpp"
#include "base/prot_mod_util.hpp"
#include "base/mass_constant.hpp"
#include "zeroptmfilter/mass_match_factory.hpp"

namespace prot {


inline std::vector<int> getScaledSrmMasses(ProteoformPtr proteo_ptr, 
                                           std::vector<double> &n_ace_shifts,
                                           double scale) {
  std::vector<int> masses = proteo_ptr->getBpSpecPtr()->getScaledPrmMasses(scale);
  for (size_t i = 0; i < n_ace_shifts.size(); i++) {
      int ace_mass = (int)std::round(- n_ace_shifts[i] * scale);
      //LOG_DEBUG("ace shift " << n_ace_shifts[i] <<  " ace mass " << ace_mass);
      masses.push_back(ace_mass);
  }
  std::sort(masses.begin(), masses.end(),std::less<int>()); 
  /*
  for (size_t i = 0; i < masses.size(); i++) {
    LOG_DEBUG("integer mass " << i << " " << masses[i]);
  }
  */
  std::vector<int> rev_masses;
  int len = masses.size();
  for (int i = len -1 ; i >= 0; i--) {
    rev_masses.push_back(masses[len-1] - masses[i]);
  }
  return rev_masses;
}

inline MassMatchPtr getMassMatchPtr(const ProteoformPtrVec &proteo_ptrs,
                                    std::vector<std::vector<int>> &mass_2d,
                                    double max_proteoform_mass, double scale, bool prm) {
  std::vector<std::vector<double>> real_shift_2d; 
  std::vector<std::vector<int>> pos_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    std::vector<double> masses;
    if (prm) {
      masses = proteo_ptrs[i]->getBpSpecPtr()->getPrmMasses();
    }
    else {
      masses = proteo_ptrs[i]->getBpSpecPtr()->getSrmMasses();
    }
    std::vector<double> shifts;
    std::vector<int> positions;
    for (size_t j = 0; j < masses.size() -1; j++) {
      shifts.push_back(-masses[j]);
      positions.push_back(j);
    }
    real_shift_2d.push_back(shifts);
    pos_2d.push_back(positions);
  }

  // for test 
  /*
  std::vector<std::vector<double>> real_shift_2d; 
  std::vector<std::vector<int>> pos_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    std::vector<double> double_masses = proteo_ptrs[i]->getBpSpecPtr()->getPrmMasses();
    std::vector<int> int_masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale);
    std::vector<double> masses;
    int len = int_masses.size();
    for (int j = 0; j < len; j++) {
      if (prm) {
        masses.push_back(int_masses[j]/scale);
      }
      else {
        masses.push_back((int_masses[len] - int_masses[j])/scale);
        std::sort(masses.begin(), masses.end(), std::less<double>());
      }
    }

    std::vector<double> shifts;
    std::vector<int> positions;
    if (prm) {
      for (size_t j = 0; j < double_masses.size() -1; j++) {
        shifts.push_back(-double_masses[j]);
        positions.push_back(j);
      }
    }
    else {
      std::sort(double_masses.begin(), double_masses.end(), std::greater<double>());
      for (size_t j = 0; j < double_masses.size() -1; j++) {
        shifts.push_back(double_masses[j] - double_masses[0]);
        positions.push_back(j);
      }
    }
    real_shift_2d.push_back(shifts);
    pos_2d.push_back(positions);
  }
  */
  
  LOG_DEBUG("pos 2d ver 1 complete");
  MassMatchPtr index_ptr(new MassMatch(mass_2d, real_shift_2d, pos_2d,
                                       max_proteoform_mass, scale));
  return index_ptr;
}


MassMatchPtr MassMatchFactory::getPrmDiagMassMatchPtr(const ProteoformPtrVec &proteo_ptrs,
                                                      double max_proteoform_mass, double scale) {
  std::vector<std::vector<int>> mass_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
      std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale);
      mass_2d.push_back(masses);
  }
  LOG_DEBUG("mass 2d ver 1 complete");
  bool prm = true;
  return getMassMatchPtr(proteo_ptrs, mass_2d, max_proteoform_mass, scale, prm);
}

MassMatchPtr MassMatchFactory::getSrmDiagMassMatchPtr(const ProteoformPtrVec &proteo_ptrs,
                                                      std::vector<std::vector<double>> &n_ace_shift_2d,
                                                      double max_proteoform_mass, double scale) {
  std::vector<std::vector<int>> mass_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    //std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledSrmMasses(scale);
    std::vector<int> masses = getScaledSrmMasses(proteo_ptrs[i], n_ace_shift_2d[i], scale);
    mass_2d.push_back(masses);
  }
  LOG_DEBUG("mass 2d ver 1 complete");
  bool prm = false;
  return getMassMatchPtr(proteo_ptrs, mass_2d, max_proteoform_mass, scale, prm);
}

inline void getPos2d (const ProteoformPtrVec &proteo_ptrs,
                      std::vector<std::vector<double>> &real_shift_2d, 
                      std::vector<std::vector<int>> &pos_2d) { 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    int n = real_shift_2d[i].size();
    std::vector<int> positions(n, 0);
    pos_2d.push_back(positions);
  }
  LOG_DEBUG("pos 2d complete");
}

MassMatchPtr MassMatchFactory::getPrmTermMassMatchPtr(const ProteoformPtrVec &proteo_ptrs, 
                                                      std::vector<std::vector<double>> &real_shift_2d,
                                                      double max_proteoform_mass, double scale) {
  std::vector<std::vector<int>> mass_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale);
    mass_2d.push_back(masses);
  }
  LOG_DEBUG("mass 2d complete");
  LOG_DEBUG("proteo num " << proteo_ptrs.size());
  LOG_DEBUG("shift num " << real_shift_2d.size());

  std::vector<std::vector<int>> pos_2d; 
  getPos2d(proteo_ptrs, real_shift_2d, pos_2d);
  MassMatchPtr index_ptr(new MassMatch(mass_2d, real_shift_2d, pos_2d,
                                       max_proteoform_mass, scale));
  return index_ptr;
}

MassMatchPtr MassMatchFactory::getSrmTermMassMatchPtr(const ProteoformPtrVec &proteo_ptrs, 
                                                      std::vector<std::vector<double>> &real_shift_2d,
                                                      std::vector<std::vector<double>> &n_ace_shift_2d,
                                                      double max_proteoform_mass, double scale) {
  std::vector<std::vector<int>> mass_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    //std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledSrmMasses(scale);
    std::vector<int> masses = getScaledSrmMasses(proteo_ptrs[i], n_ace_shift_2d[i], scale);
    mass_2d.push_back(masses);
  }
  LOG_DEBUG("mass 2d complete");
  LOG_DEBUG("proteo num " << proteo_ptrs.size());
  LOG_DEBUG("shift num " << real_shift_2d.size());

  std::vector<std::vector<int>> pos_2d; 
  getPos2d(proteo_ptrs, real_shift_2d, pos_2d);
  MassMatchPtr index_ptr(new MassMatch(mass_2d, real_shift_2d, pos_2d,
                                       max_proteoform_mass, scale));
  return index_ptr;
}

}
