#include <limits>

#include "base/logger.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

PrsmStr::PrsmStr(const std::vector<std::string> &str_vec) {
  str_vec_ = str_vec;
  std::string line = PrsmUtil::getXmlLine(str_vec_, "<spectrum_id>");
  spectrum_id_ = std::stoi(PrsmUtil::getValueStr(line));
  line = PrsmUtil::getXmlLine(str_vec_, "<seq_name>");
  seq_name_ = PrsmUtil::getValueStr(line);
  line = PrsmUtil::getXmlLine(str_vec_, "<match_fragment_num>");
  match_frag_num_ = std::stod(PrsmUtil::getValueStr(line));
  line = PrsmUtil::getXmlLine(str_vec_, "<e_value>"); 
  if (line == "") { 
    e_value_ = 0.0;
  }
  else {
    e_value_ = std::stod(PrsmUtil::getValueStr(line));
  }
  line = PrsmUtil::getXmlLine(str_vec_, "<fdr>"); 
  fdr_ = std::stod(PrsmUtil::getValueStr(line));
  //LOG_DEBUG("spectrum id " << spectrum_id_ << " match num " << match_frag_num_);
}

int getXmlLineIndex(const std::vector<std::string> &str_vec,
                    const std::string &property) {
  for (size_t i = 0; i < str_vec.size(); i++) {
    size_t found = str_vec[i].find(property);
    if (found != std::string::npos) {
      return i;
    }
  }
  return -1;
}

void PrsmStr::setFdr(double fdr) {
  int i = getXmlLineIndex(str_vec_, "fdr");
  str_vec_[i] = "<fdr>" + std::to_string(fdr) + "</fdr>";
}

void PrsmStr::setId(int id) {
  int i = getXmlLineIndex(str_vec_, "prsm_id");
  str_vec_[i] = "<prsm_id>" + std::to_string(id) + "</prsm_id>";
}

}

