#include <cstddef>
#include <limits>

#include "feature/envelope_util.hpp" 

namespace prot {

int EnvelopeUtil::getMaxPos(std::vector<double> &values) {
  int max_pos = -1;
  double max_value = - std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < values.size(); i++) {
    if (values[i] > max_value) {
      max_value = values[i];
      max_pos = i;
    }
  }
  return max_pos;
}

}

