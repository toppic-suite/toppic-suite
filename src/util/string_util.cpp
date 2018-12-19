//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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
#include <string>
#include <vector>
#include <functional>
#include <algorithm>

#include "boost/algorithm/string.hpp"

#include "util/logger.hpp"
#include "util/string_util.hpp"

namespace toppic {

namespace string_util {

std::string trim(const std::string &ori_s) {
  std::string s = ori_s;
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                  std::not1(std::ptr_fun<int, int>(std::isspace))));
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
          s.end());
  return s;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::stringstream ss(s);
  std::string item;
  std::vector<std::string> elems;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::string convertToString(double value) {
  std::stringstream stream;

  if (value < 1 && value > -1 && value !=0) {
    stream << std::scientific << std::setprecision(10);
  } else {
    stream << std::fixed << std::setprecision(10);
  }
  stream << value;
  return stream.str();
}

std::string convertToString(double value, int number) {
  std::stringstream stream;
  if (value == 0) {
    stream << std::fixed << std::setprecision(0);
  } else if (value < 0.01 && value > -0.01 && value != 0) {
    if (number > 2) {
      stream << std::scientific << std::setprecision(2);
    } else {
      stream << std::scientific << std::setprecision(number);
    }
  } else {
    stream << std::fixed << std::setprecision(number);
  }
  stream << value;
  return stream.str();
}

std::string convertToScientificStr(double value, int number) {
  std::stringstream stream;
  if (value == 0) {
    stream << std::fixed << std::setprecision(0);
  } else {
    stream << std::scientific << std::setprecision(std::min(2, number));
  }
  stream << value;
  return stream.str();
}

std::string convertToString(int value) {
  std::stringstream stream;
  stream << value;
  return stream.str();
}

std::string convertToString(bool value) {
  std::stringstream stream;
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

double convertScientificToDouble(std::string str) {
  std::stringstream ss(str);
  double d = 0;
  ss >> d;
  if (ss.fail()) {
    std::string s = "Unable to format ";
    s += str;
    s += " as a number!";
    throw (s);
  }

  return (d);
}

bool endsWith(const std::string &str, const std::string &suffix) {
  return str.size() >= suffix.size() &&
      str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

} // namespace string_util

}  // namespace toppic

