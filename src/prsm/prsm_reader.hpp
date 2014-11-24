#ifndef PROT_PRSM_READER_HPP_
#define PROT_PRSM_READER_HPP_

#include <iostream>
#include <fstream>

#include "base/xml_dom_document.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

class PrsmReader {
 public:
  PrsmReader(const std::string &file_name);

  std::vector<std::string> readOnePrsmLines();

  PrsmStrPtr readOnePrsmStr();

  PrsmPtr readOnePrsm(faidx_t *fai, const ResiduePtrVec &residue_ptr_vec);

  void close();

 private:
  std::ifstream input_;
};

typedef std::shared_ptr<PrsmReader> PrsmReaderPtr;
typedef std::vector<PrsmReaderPtr> PrsmReaderPtrVec;

} /* namespace prot */

#endif /* PROT_PRSM_READER_HPP_ */
