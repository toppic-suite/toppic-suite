#ifndef PROT_FEATURE_DETECT_MNG_HPP_
#define PROT_FEATURE_DETECT_MNG_HPP_

namespace prot {

class FeatureDetectMng {
 public:
  FeatureDetectMng() {};

  // error tolerance
  double ppo_ = 0.000015;

  int intv_width_ = 500;
};

typedef std::shared_ptr<FeatureDetectMng> FeatureDetectMngPtr;

} /* namespace */

#endif 
