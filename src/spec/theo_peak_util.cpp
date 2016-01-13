#include "spec/theo_peak_util.hpp"

namespace prot {

std::vector<double> TheoPeakUtil::getTheoMassVec (const TheoPeakPtrVec &theo_peak_list) {
  std::vector<double> masses;
  for (size_t i = 0; i < theo_peak_list.size(); i++) {
    masses.push_back(theo_peak_list[i]->getModMass());
  }
  return masses;
}

}
