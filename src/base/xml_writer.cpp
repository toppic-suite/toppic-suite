#include "base/xml_writer.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

XmlWriter::XmlWriter(const std::string &file_name,
                     const std::string &root){
  file_.open(file_name.c_str());
  root_ = root;
  LOG_DEBUG("file_name " << file_name);
  file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  if(root_.compare("")!=0){
    file_ << ("<"+root_+">");
  }
  XmlDOMImpl* impl = XmlDOMImplFactory::getXmlDOMImplInstance();
  doc_ = new XmlDOMDocument(impl->createDoc(root_.compare("")!=0?root_:"ROOT"));
  serializer_ = impl->createSerializer();
}

XmlWriter::~XmlWriter(){
  serializer_->release();
  delete doc_;
}

void XmlWriter::write(xercesc::DOMElement* element){
  std::string str = XmlDomUtil::writeToString(serializer_, element);
  XmlDomUtil::writeToStreamByRemovingDoubleLF(file_, str);
  element->release();
}

void XmlWriter::close(){
  if(root_.compare("")!=0){
    file_ << "</"+root_+">" << std::endl;
  }
  file_.close();
}

} /* namespace prot */
