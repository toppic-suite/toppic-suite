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

class PrsmWriter {
public:
 PrsmWriter(std::string file_name);
 ~PrsmWriter();
 void close();
 void write(PrsmPtr prsm_ptr);
 void writeVector(PrsmPtrVec &prsms);
 void writeVector2D(PrsmPtrVec2D &prsms);
 void writeVector3D(PrsmPtrVec3D &prsms);
private:
  //XmlDOMDocument* doc_;
  std::ofstream file_;
};

typedef std::shared_ptr<PrsmWriter> PrsmWriterPtr;
typedef std::vector<PrsmWriterPtr> PrsmWriterPtrVec;

} /* namespace prot */

#endif /* PRSM_WRITER_HPP_ */
