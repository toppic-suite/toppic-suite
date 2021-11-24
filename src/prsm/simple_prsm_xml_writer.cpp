//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include "common/util/file_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/xml/xml_dom_impl.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"

namespace toppic {

SimplePrsmXmlWriter::SimplePrsmXmlWriter(const std::string &file_name) {
  file_.open(file_name.c_str());
  file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  file_ << "<simple_prsm_list>" << std::endl;
  XmlDOMImpl* impl = XmlDOMImplFactory::getXmlDOMImplInstance();
  doc_ = new XmlDOMDocument(impl->createDoc("simple_prsm_list"));
  serializer_ = impl->createSerializer();

  file_name_ = file_util::basename(file_name) + ".msalign";
}

SimplePrsmXmlWriter::~SimplePrsmXmlWriter() {
  serializer_->release();
  delete doc_;
  if (file_.is_open()) {
    close();
  }
}

void SimplePrsmXmlWriter::close() {
  file_ << "</simple_prsm_list>" << std::endl;
  file_.close();
}

void SimplePrsmXmlWriter::write(SimplePrsmStrPtr prsm_str_ptr) {
  std::vector<std::string> strs = prsm_str_ptr->getStrVec();
  for (size_t i = 0; i < strs.size(); i++) {
    file_ << strs[i] << std::endl;
  }
}

void SimplePrsmXmlWriter::write(const SimplePrsmPtrVec &simple_prsm_ptrs) {
  for (size_t i = 0; i < simple_prsm_ptrs.size(); i++) {
    write(simple_prsm_ptrs[i]);
  }
}

void SimplePrsmXmlWriter::write(SimplePrsmPtr simple_prsm_ptr) {
  XmlDOMElement * element = simple_prsm_ptr->toXml(doc_);
  std::string str = xml_dom_util::writeToString(serializer_, element);
  xml_dom_util::writeToStreamByRemovingDoubleLF(file_, str);
  element->release();
}

}  // namespace toppic
