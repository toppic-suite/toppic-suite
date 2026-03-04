//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string.hpp>

#include "common/util/str_util.hpp"

namespace toppic {

namespace str_util {

void trim(std::string &s) {
  boost::algorithm::trim(s);
}

std::vector<std::string> split(const std::string &s, const std::string &delim) {
  std::vector<std::string> strs;
  boost::split(strs, s, boost::is_any_of(delim));
  return strs;
}

std::string toString(bool value) {
  return value ? "true" : "false";
}

std::string toString(int value) {
  return std::to_string(value);
}

std::string toString(size_t value) {
  return std::to_string(value);
}

std::string toString(double value) {
  std::stringstream stream;
  if (value != 0 && value < 1 && value > -1) {
    stream << std::scientific << std::setprecision(10);
  } else {
    stream << std::fixed << std::setprecision(10);
  }
  stream << value;
  return stream.str();
}

std::string evalueToString(double value, int precision) {
  std::stringstream stream;
  if (value == 0) {
    stream << std::fixed << std::setprecision(0);
  } else if (value < 0.01 && value > -0.01) {
    if (precision > 2) {
      stream << std::scientific << std::setprecision(2);
    } else {
      stream << std::scientific << std::setprecision(precision);
    }
  } else {
    stream << std::fixed << std::setprecision(precision);
  }
  stream << value;
  return stream.str();
}

std::string confToString(double value, int precision) {
  return evalueToString(value, precision);
}

std::string fixedToString(double value, int precision) {
  std::stringstream stream;
  if (value == 0) {
    stream << std::fixed << std::setprecision(0);
  } else {
    stream << std::fixed << std::setprecision(precision);
  }
  stream << value;
  return stream.str();
}

std::string toScientificStr(double value, int precision) {
  std::stringstream stream;
  if (value == 0) {
    stream << std::fixed << std::setprecision(0);
  } else {
    stream << std::scientific << std::setprecision(precision);
  }
  stream << value;
  return stream.str();
}

std::string rmComment(const std::string &ori_s, const std::string &comment) {
  std::string s = ori_s;
  std::string::size_type i = s.find(comment);
  if (i != std::string::npos) s.erase(i);
  boost::trim_right(s);
  return s;
}

double scientificToDouble(const std::string &str) {
  std::stringstream ss(str);
  double d = 0;
  ss >> d;
  if (ss.fail()) {
    throw std::invalid_argument("Cannot convert \"" + str + "\" to double!");
  }
  return d;
}

bool endsWith(const std::string &str, const std::string &suffix) {
  return str.size() >= suffix.size() &&
      str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

}  // namespace str_util

}  // namespace toppic
