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
