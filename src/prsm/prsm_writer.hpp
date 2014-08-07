#ifndef PRSM_WRITER_HPP_
#define PRSM_WRITER_HPP_

#include <iostream>
#include <fstream>

#include "base/xml_dom_document.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class PrsmWriter {
public:
 PrsmWriter(const std::string &file_name);

 void close();
 void write(PrsmPtr prsm_ptr);
 void writeVector(const PrsmPtrVec &prsm_ptrs);
 void writeVector2D(const PrsmPtrVec2D &prsm_ptrs);
 void writeVector3D(const PrsmPtrVec3D &prsm_ptrs);
private:
  //XmlDOMDocument* doc_;
  std::ofstream file_;
};

typedef std::shared_ptr<PrsmWriter> PrsmWriterPtr;
typedef std::vector<PrsmWriterPtr> PrsmWriterPtrVec;

} /* namespace prot */

#endif /* PRSM_WRITER_HPP_ */
