#include <iostream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/prot_mod.hpp"
#include "base/prot_mod_util.hpp"
#include "base/mass_constant.hpp"
#include "zeroptmfilter/mass_match_factory.hpp"

namespace prot {

MassMatchPtr MassMatchFactory::getMassMatchPtr(const ProteoformPtrVec &proteo_ptrs, 
                                               double max_proteoform_mass, double scale, bool rev) {
  std::vector<std::vector<int>> mass_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    if (!rev) {
      std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale);
      mass_2d.push_back(masses);
    }
    else {
      std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledSrmMasses(scale);
      mass_2d.push_back(masses);
    }
  }
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
      std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledSrmMasses(scale);
      mass_2d.push_back(masses);
    }
  }
  std::vector<std::vector<int>> pos_2d; 
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    std::vector<int> positions(real_shift_2d[i].size(), 0);
    pos_2d.push_back(positions);
  }
  MassMatchPtr index_ptr(new MassMatch(mass_2d, real_shift_2d, pos_2d,
                                       max_proteoform_mass, scale));
  return index_ptr;
}

}
