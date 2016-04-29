#include <iostream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/prot_mod.hpp"
#include "base/prot_mod_util.hpp"
#include "base/mass_constant.hpp"
#include "zeroptmfilter/mass_match_factory.hpp"

namespace prot {


std::vector<int> getScaledSrmMasses(ProteoformPtr proteo_ptr, double scale) {
  std::vector<int> masses = proteo_ptr->getBpSpecPtr()->getScaledPrmMasses(scale);
  std::vector<int> rev_masses;
  int len = masses.size();
  for (int i = len -1 ; i >= 0; i--) {
    rev_masses.push_back(masses[len-1] - masses[i]);
  }
  return rev_masses;
}


MassMatchPtr MassMatchFactory::getMassMatchPtr(const ProteoformPtrVec &proteo_ptrs, 
                                               double max_proteoform_mass, double scale, bool rev) {
  std::vector<std::vector<int>> mass_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    if (!rev) {
      std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale);
      mass_2d.push_back(masses);
    }
    else {
      //std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledSrmMasses(scale);
      std::vector<int> masses = getScaledSrmMasses(proteo_ptrs[i], scale);
      mass_2d.push_back(masses);
    }
  }
  LOG_DEBUG("mass 2d ver 1 complete");
  std::vector<std::vector<double>> real_shift_2d; 
  std::vector<std::vector<int>> pos_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    std::vector<double> masses;
    if (!rev) {
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
  LOG_DEBUG("pos 2d ver 1 complete");
  MassMatchPtr index_ptr(new MassMatch(mass_2d, real_shift_2d, pos_2d,
                                       max_proteoform_mass, scale));
  return index_ptr;
}

MassMatchPtr MassMatchFactory::getMassMatchPtr(const ProteoformPtrVec &proteo_ptrs, 
                                               std::vector<std::vector<double>> &real_shift_2d,
                                               double max_proteoform_mass, double scale, bool rev) {
  std::vector<std::vector<int>> mass_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    if (!rev) {
      std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale);
      mass_2d.push_back(masses);
    }
    else {
      //std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledSrmMasses(scale);
      std::vector<int> masses = getScaledSrmMasses(proteo_ptrs[i], scale);
      mass_2d.push_back(masses);
    }
  }
  LOG_DEBUG("mass 2d complete");
  LOG_DEBUG("proteo num " << proteo_ptrs.size());
  LOG_DEBUG("shift num " << real_shift_2d.size());

  std::vector<std::vector<int>> pos_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    int n = real_shift_2d[i].size();
    std::vector<int> positions(n, 0);
    pos_2d.push_back(positions);
  }
  LOG_DEBUG("pos 2d complete");
  MassMatchPtr index_ptr(new MassMatch(mass_2d, real_shift_2d, pos_2d,
                                       max_proteoform_mass, scale));
  return index_ptr;
}

}
