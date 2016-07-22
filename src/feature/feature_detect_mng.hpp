#ifndef PROT_FEATURE_DETECT_MNG_HPP_
#define PROT_FEATURE_DETECT_MNG_HPP_

#include <vector>
#include "spec/peak_tolerance.hpp"

namespace prot {

class FeatureDetectMng {
 public:
  FeatureDetectMng();

  std::vector<double> getExtMasses(double mass);

  PeakTolerancePtr peak_tolerance_ptr_;

  std::vector<double> ext_offsets_;

  double extend_min_mass_ = 5000;

  int intv_width_ = 500;

};

typedef std::shared_ptr<FeatureDetectMng> FeatureDetectMngPtr;

} /* namespace */

#endif 
