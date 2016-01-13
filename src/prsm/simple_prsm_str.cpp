#include "base/logger.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/simple_prsm_str.hpp"

namespace prot {

SimplePrsmStr::SimplePrsmStr(const std::vector<std::string> &str_vec) {
  str_vec_ = str_vec;
  std::string line = PrsmUtil::getXmlLine(str_vec_,"spectrum_id"); 
  spectrum_id_ = std::stoi(PrsmUtil::getValueStr(line));
  line = PrsmUtil::getXmlLine(str_vec_,"score"); 
  score_ = std::stod(PrsmUtil::getValueStr(line));
  //LOG_DEBUG("spectrum id " << spectrum_id_ << " score " << score_);
}

}
