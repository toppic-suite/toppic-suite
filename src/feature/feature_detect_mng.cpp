#include "base/mass_constant.hpp"
#include "feature/feature_detect_mng.hpp"

namespace prot {

FeatureDetectMng::FeatureDetectMng() {
  double ppo = 0.000015;
  bool use_min_tolerance = true;
  double min_tolerance = 0.01;
  peak_tolerance_ptr_ = PeakTolerancePtr(
      new PeakTolerance(ppo, use_min_tolerance, min_tolerance));

  // extend sp parameter 
  double IM = MassConstant::getIsotopeMass();
  // the set of offsets used to expand the monoisotopic mass list 
  std::vector<double> offsets {{0, -IM, IM, -2 * IM, 2 * IM}};
  ext_offsets_ = offsets;
}

std::vector<double> FeatureDetectMng::getExtMasses(double mass) {
  std::vector<double> result;
  if (mass < extend_min_mass_) {
    result.push_back(mass);
  }
  else {
    for (size_t i = 0; i < ext_offsets_.size(); i++) {
      double new_mass = mass + ext_offsets_[i];
      result.push_back(new_mass);
    }
  }
  return result;
}

} /* namespace */

