#ifndef PROT_PRSM_PRSM_XML_WRITER_HPP_
#define PROT_PRSM_PRSM_XML_WRITER_HPP_

#include <iostream>
#include <fstream>

#include "base/xml_dom_document.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

class PrsmXmlWriter {
public:
 PrsmXmlWriter(const std::string &file_name);

 void close();

 void write(PrsmStrPtr prsm_str_ptr);
 void writeVector(const PrsmStrPtrVec &prsm_str_ptr_vec);

 void write(PrsmPtr prsm_ptr);
 void writeVector(const PrsmPtrVec &prsm_ptrs);
 void writeVector2D(const PrsmPtrVec2D &prsm_ptrs);
 void writeVector3D(const PrsmPtrVec3D &prsm_ptrs);
private:
  //XmlDOMDocument* doc_;
  std::ofstream file_;
};

typedef std::shared_ptr<PrsmXmlWriter> PrsmXmlWriterPtr;
typedef std::vector<PrsmXmlWriterPtr> PrsmXmlWriterPtrVec;

} /* namespace prot */

#endif /* PRSM_WRITER_HPP_ */
