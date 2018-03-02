//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#include <cmath>
#include <utility>
#include <algorithm>
#include <vector>

#include "base/mass_constant.hpp"
#include "spec/extend_ms.hpp"

namespace prot {

namespace extend_ms {

std::vector<double> getExtendMassVec(ExtendMsPtr extend_ms_ptr) {
  std::vector<double> masses;
  ExtendPeakPtrVec peak_ptr_list = extend_ms_ptr->getPeakPtrVec();
  for (size_t i = 0; i < peak_ptr_list.size(); i++) {
    masses.push_back(peak_ptr_list[i]->getPosition());
  }
  return masses;
}

inline bool massErrorUp(const std::pair<int, int> &a, const std::pair<int, int> b) {
  return a.first < b.first;
}

std::vector<std::pair<int, int>> getExtendIntMassErrorList(const ExtendMsPtrVec &ext_ms_ptr_vec,
                                                           bool pref, double scale) {
  std::vector<std::pair<int, int>> mass_errors;
  for (size_t i = 0; i < ext_ms_ptr_vec.size(); i++) {
    ExtendMsPtr ext_ms_ptr = ext_ms_ptr_vec[i];
    double shift = 0;
    if (pref) {
      shift = ext_ms_ptr->getMsHeaderPtr()->getActivationPtr()->getNShift();
    } else {
      shift = ext_ms_ptr->getMsHeaderPtr()->getActivationPtr()->getCShift()
          + mass_constant::getWaterMass();
    }

    std::pair<int, int> last_mass_error(-1, 0);
    for (size_t j = 0; j < ext_ms_ptr->size(); j++) {
      int m = static_cast<int>(std::round((ext_ms_ptr->getPeakPtr(j)->getPosition() - shift) * scale));
      int e = static_cast<int>(std::ceil(ext_ms_ptr->getPeakPtr(j)->getOrigTolerance() * scale));
      std::pair<int, int> cur_m_e(m, e);
      if (cur_m_e.first != last_mass_error.first) {
        mass_errors.push_back(cur_m_e);
        last_mass_error = cur_m_e;
      } else if (cur_m_e.second > last_mass_error.second) {
        mass_errors.pop_back();
        mass_errors.push_back(cur_m_e);
        last_mass_error = cur_m_e;
      }
    }
  }

  std::sort(mass_errors.begin(), mass_errors.end(), massErrorUp);
  return mass_errors;
}

}  // namespace extend_ms

}  // namespace prot
