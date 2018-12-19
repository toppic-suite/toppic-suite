//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#include "xml/xml_writer.hpp"
#include "xml/xml_dom_util.hpp"

namespace toppic {

XmlWriter::XmlWriter(const std::string &file_name, const std::string &root) {
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

XmlWriter::~XmlWriter() {
  serializer_->release();
  delete doc_;
}

void XmlWriter::write(xercesc::DOMElement* element) {
  std::string str = xml_dom_util::writeToString(serializer_, element);
  xml_dom_util::writeToStreamByRemovingDoubleLF(file_, str);
  element->release();
}

void XmlWriter::write_str(const std::string & str) {
  file_ << str << std::endl;
}

void XmlWriter::close() {
  if(root_.compare("")!=0){
    file_ << "</"+root_+">" << std::endl;
  }
  file_.close();
}

} /* namespace toppic */
