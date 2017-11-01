//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include "base/logger.hpp"
#include "base/string_util.hpp"

namespace prot {

std::string StringUtil::trim(const std::string &ori_s) {
  std::string s = ori_s;
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                  std::not1(std::ptr_fun<int, int>(std::isspace))));
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
          s.end());
  return s;
}

std::vector<std::string> StringUtil::split(const std::string &s, char delim) {
  std::stringstream ss(s);
  std::string item;
  std::vector<std::string> elems;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::string StringUtil::convertToString(double value) {
  std::stringstream stream;

  if (value < 1 && value > -1 && value !=0) {
    stream << std::scientific << std::setprecision(10);
  } else {
    stream << std::fixed << std::setprecision(10);
  }
  stream << value;
  return stream.str();
}

std::string StringUtil::convertToString(double value, int number) {
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

std::string StringUtil::convertToScientificStr(double value, int number) {
  std::stringstream stream;
  if (value == 0) {
    stream << std::fixed << std::setprecision(0);
  } else {
    stream << std::scientific << std::setprecision(std::min(2, number));
  }
  stream << value;
  return stream.str();
}

std::string StringUtil::convertToString(int value) {
  std::stringstream stream;
  stream << value;
  return stream.str();
}

std::string StringUtil::convertToString(bool value) {
  std::stringstream stream;
  stream << value;
  return stream.str();
}

std::string StringUtil::rmComment(const std::string &ori_s, const std::string &comment) {
  std::string s = ori_s;
  std::string::size_type i = s.find(comment);

  if (i != std::string::npos) s.erase(i);

  boost::trim_right(s);

  return s;
}

}  // namespace prot

