#ifndef PROT_BASE_STRING_UTIL_HPP_
#define PROT_BASE_STRING_UTIL_HPP_

#include <sstream>
#include <string>
#include <vector>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

#include "base/logger.hpp"

namespace prot {

class StringUtil {
 public:
  static std::string trim(std::string &ori_s);

  static std::vector<std::string> split(const std::string &ori_s, char delim);

  static std::string convertToString(double value);

  static std::string convertToString(double value, int number);

  static std::string convertToString(int value);

  static std::string convertToString(bool value);
};

}

#endif
