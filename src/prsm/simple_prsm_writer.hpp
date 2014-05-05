/*
 * simple_prsm_writer.hpp
 *
 *  Created on: Dec 9, 2013
 *      Author: xunlikun
 */

#ifndef PROT_SIMPLE_PRSM_WRITER_HPP_
#define PROT_SIMPLE_PRSM_WRITER_HPP_

#include "prsm/simple_prsm.hpp"

namespace prot {

class SimplePrsmWriter {
 public:
  SimplePrsmWriter();
  SimplePrsmWriter(std::string file_name);
  ~SimplePrsmWriter();
  void close();
  void write(SimplePrsmPtrVec simple_prsms);

 private:
  xercesc::DOMLSSerializer* serializer_;
  XmlDOMDocument* doc_;
  std::ofstream file_;
};

} /* namespace prot */

#endif /* SIMPLE_PRSM_WRITER_HPP_ */
