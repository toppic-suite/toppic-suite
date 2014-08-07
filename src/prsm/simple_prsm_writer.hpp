#ifndef PROT_SIMPLE_PRSM_WRITER_HPP_
#define PROT_SIMPLE_PRSM_WRITER_HPP_

#include "prsm/simple_prsm.hpp"

namespace prot {

class SimplePrsmWriter {
 public:
  SimplePrsmWriter();
  SimplePrsmWriter(const std::string &file_name);
  ~SimplePrsmWriter();
  void close();
  void write(const SimplePrsmPtrVec &simple_prsm_ptrs);

 private:
  xercesc::DOMLSSerializer* serializer_;
  XmlDOMDocument* doc_;
  std::ofstream file_;
};

} /* namespace prot */

#endif /* SIMPLE_PRSM_WRITER_HPP_ */
