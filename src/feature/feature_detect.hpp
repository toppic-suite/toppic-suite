#ifndef PROT_FEATURE_DETECTION_HPP_
#define PROT_FEATURE_DETECTION_HPP_

#include <string>

namespace prot {

class FeatureDetect {
 public:
  static void process(std::string &input_file_name);
};

} /* namespace_prot */

#endif
