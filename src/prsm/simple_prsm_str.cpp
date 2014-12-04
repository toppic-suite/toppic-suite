#include "base/logger.hpp"
#include "prsm/simple_prsm_str.hpp"

namespace prot {

std::string getValueStr(std::string line) {
  int start = line.find(">");
  int end = line.find("<", start);
  std::string num_str = line.substr(start + 1, end - start - 1);
  //LOG_DEBUG(line << "  Num str: " << num_str);
  return num_str;
}

std::string getXmlLine(const std::vector<std::string> &str_vec,
                       const std::string &property) {
  for (size_t i = 0; i < str_vec.size(); i++) {
    size_t found = str_vec[i].find(property);
    if (found != std::string::npos) {
      return str_vec[i];
    }
  }
  return "";
}

SimplePrsmStr::SimplePrsmStr(const std::vector<std::string> &str_vec) {
  str_vec_ = str_vec;
  std::string line = getXmlLine(str_vec_,"spectrum_id"); 
  spectrum_id_ = std::stoi(getValueStr(line));
  line = getXmlLine(str_vec_,"score"); 
  score_ = std::stod(getValueStr(line));
  //LOG_DEBUG("spectrum id " << spectrum_id_ << " score " << score_);
}


}

