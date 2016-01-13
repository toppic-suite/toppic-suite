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

