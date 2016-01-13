#include <iostream>
#include <algorithm>

#include "base/xml_dom_document.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_util.hpp"
#include "base/file_util.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"

namespace prot {

SimplePrsmXmlWriter::SimplePrsmXmlWriter(const std::string &file_name){
    file_.open(file_name.c_str());
    file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    file_ << "<simple_prsm_list>" << std::endl;
    XmlDOMImpl* impl = XmlDOMImplFactory::getXmlDOMImplInstance();
    doc_ = new XmlDOMDocument(impl->createDoc("simple_prsm_list"));
    serializer_ = impl->createSerializer();
}

SimplePrsmXmlWriter::~SimplePrsmXmlWriter(){
    serializer_->release();
    delete doc_;
}

void SimplePrsmXmlWriter::close(){
  file_ << "</simple_prsm_list>" << std::endl;
  file_.close();
}

void SimplePrsmXmlWriter::write(SimplePrsmStrPtr prsm_str_ptr) {
  std::vector<std::string> strs = prsm_str_ptr->getStrVec();
  for(size_t i = 0; i < strs.size(); i++) {
    file_ << strs[i] << std::endl;
  }
}

void SimplePrsmXmlWriter::write(const SimplePrsmPtrVec &simple_prsm_ptrs){
  for(size_t i=0;i<simple_prsm_ptrs.size();i++){
    xercesc::DOMElement* element = simple_prsm_ptrs[i]->toXml(doc_);
    std::string str = XmlDomUtil::writeToString(serializer_, element);
    XmlDomUtil::writeToStreamByRemovingDoubleLF(file_, str);
    element->release();
  }
}

} /* namespace prot */
