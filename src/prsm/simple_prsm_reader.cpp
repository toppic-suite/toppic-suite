#include <xercesc/framework/MemBufInputSource.hpp>

#include "base/logger.hpp"
#include "base/string_util.hpp"
#include "prsm/simple_prsm_reader.hpp"

namespace prot {

SimplePrsmReader::SimplePrsmReader(const std::string &file_name) {
  input_.open(file_name.c_str(), std::ios::in);
}

std::vector<std::string> SimplePrsmReader::readOnePrsmLines() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    line = trim(line);
    if (line ==  "<simple_prsm>") {
      line_list.push_back(line);
    }
    else if (line == "</simple_prsm>") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    }
    else if (line == "") {
      continue;
    }
    else {
      if (line_list.size() > 0) {
        line_list.push_back(line);
      }
    }
  }
  return line_list;
}

SimplePrsmStrPtr SimplePrsmReader::readOnePrsmStr() {
  std::vector<std::string> prsm_str_vec = readOnePrsmLines();
  if (prsm_str_vec.size() == 0) {
    return SimplePrsmStrPtr(nullptr);
  }
  return SimplePrsmStrPtr(new SimplePrsmStr(prsm_str_vec));
}


void SimplePrsmReader::close() {
  input_.close();
}


} /* namespace prot */
