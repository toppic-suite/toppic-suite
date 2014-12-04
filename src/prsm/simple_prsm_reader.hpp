#ifndef PROT_SIMPLE_PRSM_READER_HPP_
#define PROT_SIMPLE_PRSM_READER_HPP_

#include <iostream>
#include <fstream>

#include "base/xml_dom_document.hpp"
#include "prsm/prsm.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_str.hpp"

namespace prot {

class SimplePrsmReader {
 public:
  SimplePrsmReader(const std::string &file_name);

  std::vector<std::string> readOnePrsmLines();

  SimplePrsmStrPtr readOnePrsmStr();

  SimplePrsmPtr readOnePrsm();

  void close();

 private:
  std::ifstream input_;
};

typedef std::shared_ptr<SimplePrsmReader> SimplePrsmReaderPtr;
typedef std::vector<SimplePrsmReaderPtr> SimplePrsmReaderPtrVec;

} /* namespace prot */

#endif /* PROT_SIMPLE_PRSM_READER_HPP_ */
