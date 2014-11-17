#include "base/logger.hpp"
#include "prsm/simple_prsm_str.hpp"

namespace prot {

std::string getNumStr(std::string line) {
  int start = line.find(">");
  int end = line.find("<", start);
  std::string num_str = line.substr(start + 1, end - start - 1);
  //LOG_DEBUG(line << "  Num str: " << num_str);
  return num_str;
}

SimplePrsmStr::SimplePrsmStr(const std::vector<std::string> &str_vec) {
  str_vec_ = str_vec;
  spectrum_id_ = std::stoi(getNumStr(str_vec_[1]));
  score_ = std::stod(getNumStr(str_vec_[6]));
  //LOG_DEBUG("spectrum id " << spectrum_id_ << " score " << score_);
}


}

