// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <iomanip>

#include "base/logger.hpp"
#include "base/string_util.hpp"

namespace prot {

std::string StringUtil::trim(std::string &ori_s) {
  std::string s = ori_s;
  s.erase(s.begin(), 
          std::find_if(s.begin(), s.end(), 
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

  if(value < 1 && value > -1 && value !=0){
    stream << std::scientific << std::setprecision(10);
  }
  else{
    stream << std::fixed<<std::setprecision(10);
  }
  stream << value;
  return stream.str();
}

std::string StringUtil::convertToString(double value, int number) {
  std::stringstream stream;
  if(value ==0) {
    stream << std::fixed << std::setprecision(0);
  }
  else if (value < 0.01 && value > -0.01 && value != 0) {
    if(number>2){
      stream << std::scientific << std::setprecision(2);
    }
    else{
      stream << std::scientific << std::setprecision(number);
    }
  } else {
    stream << std::fixed << std::setprecision(number);
  }
  stream << value;
  return stream.str();
}

std::string StringUtil::convertToString(int value){
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

  if (i != std::string::npos)                                                   
    s.erase(i);

  return s;
}

}

