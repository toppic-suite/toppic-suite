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

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "common/util/str_util.hpp"

namespace toppic {

namespace str_util {

namespace {

// True for non-whitespace characters. The cast to unsigned char avoids the
// undefined behavior of passing a negative char to std::isspace.
bool isNotSpace(unsigned char c) { return std::isspace(c) == 0; }

}  // namespace

void trim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), isNotSpace));
  s.erase(std::find_if(s.rbegin(), s.rend(), isNotSpace).base(), s.end());
}

// Splits on any single character contained in delim, preserving empty fields.
// This matches the previous boost::split + boost::is_any_of behavior, e.g.
// "a,,b" -> {"a", "", "b"}, ",a" -> {"", "a"}, and "" -> {""}.
std::vector<std::string> split(const std::string &s, const std::string &delim) {
  std::vector<std::string> strs;
  std::string token;
  for (char c : s) {
    if (delim.find(c) != std::string::npos) {
      strs.push_back(token);
      token.clear();
    } else {
      token.push_back(c);
    }
  }
  strs.push_back(token);
  return strs;
}

std::string toString(bool value) {
  return value ? "true" : "false";
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
  // Right-trim trailing whitespace (replaces boost::trim_right).
  s.erase(std::find_if(s.rbegin(), s.rend(), isNotSpace).base(), s.end());
  return s;
}

}  // namespace str_util

}  // namespace toppic
