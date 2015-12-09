#include "base/logger.hpp"
#include "prsm/prsm_util.hpp"

namespace prot {

std::string PrsmUtil::getValueStr(std::string line) {
  int start = line.find(">");
  int end = line.find("<", start);
  std::string num_str = line.substr(start + 1, end - start - 1);
  //LOG_DEBUG(line << "  Num str: " << num_str);
  return num_str;
}

std::string PrsmUtil::getXmlLine(const std::vector<std::string> &str_vec,
                                       const std::string &property) {
  for (size_t i = 0; i < str_vec.size(); i++) {
    size_t found = str_vec[i].find(property);
    if (found != std::string::npos) {
      return str_vec[i];
    }
  }
  return "";
}

}

