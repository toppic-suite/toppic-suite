/*
 * prsm_writer.hpp
 *
 *  Created on: Dec 30, 2013
 *      Author: xunlikun
 */

#ifndef PRSM_WRITER_HPP_
#define PRSM_WRITER_HPP_

#include <iostream>
#include <fstream>

#include "base/xml_dom_document.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class PrSMWriter {
public:
 PrSMWriter(std::string file_name);
 ~PrSMWriter();
 void write(PrSMPtr prsm_ptr);
 void writeVector(PrSMPtrVec &prsms);

private:
  xercesc::DOMLSSerializer* serializer_;
  XmlDOMDocument* doc_;
  std::ofstream file_;
};

} /* namespace prot */

#endif /* PRSM_WRITER_HPP_ */
