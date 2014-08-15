#include "base/logger.hpp"
#include "prsm/prsm_writer.hpp"

namespace prot {

PrsmWriter::PrsmWriter(const std::string &file_name) {
  file_.open(file_name.c_str());
  LOG_DEBUG("file_name " << file_name);
  file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  file_ << "<prsm_list>";
}

void PrsmWriter::close(){
  file_ << "</prsm_list>" << std::endl;
  file_.close();
}


void PrsmWriter::write(PrsmPtr prsm_ptr) {
  if(prsm_ptr!=nullptr){
    XmlDOMImpl* impl = XmlDOMImplFactory::getXmlDOMImplInstance();
    xercesc::DOMLSSerializer* serializer = impl->createSerializer();
    XmlDOMDocument doc (impl->createDoc("prsm_list"));
    xercesc::DOMElement* element = prsm_ptr->toXmlElement(&doc);
    //LOG_DEBUG("Element generated");
    std::string str = writeToString(serializer, element);
    //LOG_DEBUG("String generated");
    writeToStreamByRemovingDoubleLF(file_, str);
    element->release();
    serializer->release();
  }
}

void PrsmWriter::writeVector(const PrsmPtrVec &prsm_ptrs) {
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    write(prsm_ptrs[i]);
  }
}

void PrsmWriter::writeVector2D(const PrsmPtrVec2D &prsm_ptrs){
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    writeVector(prsm_ptrs[i]);
  }
}

void PrsmWriter::writeVector3D(const PrsmPtrVec3D &prsm_ptrs){
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    writeVector2D(prsm_ptrs[i]);
  }
}

} /* namespace prot */
