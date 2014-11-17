
#include "base/logger.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

std::string getNumStr(const std::string &line) {
  int start = line.find(">");
  int end = line.find("<", start);
  std::string num_str = line.substr(start + 1, end - start - 1);
  //LOG_DEBUG(line << "  Num str: " << num_str);
  return num_str;
}

PrsmStr::PrsmStr(const std::vector<std::string> &str_vec) {
  str_vec_ = str_vec;
  spectrum_id_ = std::stoi(getNumStr(str_vec_[2]));
  match_frag_num_ = std::stod(getNumStr(str_vec_[10]));
  //LOG_DEBUG("spectrum id " << spectrum_id_ << " match num " << match_frag_num_);
}


}

