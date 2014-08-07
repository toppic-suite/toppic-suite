#include <iostream>
#include "prsm/simple_prsm_writer.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom.hpp"

namespace prot {

SimplePrsmWriter::SimplePrsmWriter(const std::string &file_name){
    file_.open(file_name.c_str());
    file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    file_ << "<simple_prsm_list>" << std::endl;
    XmlDOMImpl* impl = XmlDOMImplFactory::getXmlDOMImplInstance();
    doc_ = new XmlDOMDocument(impl->createDoc("simple_prsm_list"));
    serializer_ = impl->createSerializer();
}

SimplePrsmWriter::~SimplePrsmWriter(){
    serializer_->release();
    delete doc_;
}

void SimplePrsmWriter::close(){
  file_ << "</simple_prsm_list>" << std::endl;
  file_.close();
}

void SimplePrsmWriter::write(const SimplePrsmPtrVec &simple_prsm_ptrs){
  for(size_t i=0;i<simple_prsm_ptrs.size();i++){
    xercesc::DOMElement* element = simple_prsm_ptrs[i]->toXml(doc_);
    std::string str = writeToString(serializer_, element);
    writeToStreamByRemovingDoubleLF(file_, str);
    element->release();
  }
}

} /* namespace prot */
