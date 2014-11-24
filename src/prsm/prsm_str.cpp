#include <limits>

#include "base/logger.hpp"
#include "prsm/simple_prsm_str.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

PrsmStr::PrsmStr(const std::vector<std::string> &str_vec) {
  str_vec_ = str_vec;
  std::string line = getXmlLine(str_vec_, "spectrum_id");
  spectrum_id_ = std::stoi(getValueStr(line));
  line = getXmlLine(str_vec_, "db_id");
  db_id_ = std::stoi(getValueStr(line));
  line = getXmlLine(str_vec_, "match_fragment_num");
  match_frag_num_ = std::stod(getValueStr(line));
  line = getXmlLine(str_vec_, "e_value"); 
  if (line == "") { 
    e_value_ = std::numeric_limits<double>::max();
  }
  else {
    e_value_ = std::stod(getValueStr(line));
  }
  //LOG_DEBUG("spectrum id " << spectrum_id_ << " match num " << match_frag_num_);
}


}

