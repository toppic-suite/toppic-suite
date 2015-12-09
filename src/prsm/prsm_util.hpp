#ifndef PROT_PRSM_PRSM_UTIL_HPP_
#define PROT_PRSM_PRSM_UTIL_HPP_

#include <memory>
#include <vector>
#include <string>

namespace prot {

class PrsmUtil {
 public:
  static std::string getValueStr(std::string line);

  static std::string getXmlLine(const std::vector<std::string> &str_vec,
                                const std::string &property);
};

}
#endif

