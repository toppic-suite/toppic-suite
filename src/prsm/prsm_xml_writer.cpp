//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include "base/logger.hpp"
#include "base/xml_dom_util.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace prot {

PrsmXmlWriter::PrsmXmlWriter(const std::string &file_name) {
  file_.open(file_name.c_str());
  LOG_DEBUG("file_name " << file_name);
  file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  file_ << "<prsm_list>" << std::endl;
}

void PrsmXmlWriter::close(){
  file_ << "</prsm_list>" << std::endl;
  file_.close();
}

void PrsmXmlWriter::write(PrsmStrPtr prsm_str_ptr) {
  std::vector<std::string> strs = prsm_str_ptr->getStrVec();
  for(size_t i = 0; i < strs.size(); i++) {
    file_ << strs[i] << std::endl;
  }
}

void PrsmXmlWriter::writeVector(const PrsmStrPtrVec &prsm_str_ptr_vec) {
  for(size_t i = 0; i < prsm_str_ptr_vec.size(); i++) {
    write(prsm_str_ptr_vec[i]);
  }
}


void PrsmXmlWriter::write(PrsmPtr prsm_ptr) {
  if(prsm_ptr!=nullptr){
    XmlDOMImpl* impl = XmlDOMImplFactory::getXmlDOMImplInstance();
    xercesc::DOMLSSerializer* serializer = impl->createSerializer();
    XmlDOMDocument doc (impl->createDoc("prsm_list"));
    xercesc::DOMElement* element = prsm_ptr->toXmlElement(&doc);
    //LOG_DEBUG("Element generated");
    std::string str = XmlDomUtil::writeToString(serializer, element);
    //LOG_DEBUG("String generated");
    XmlDomUtil::writeToStreamByRemovingDoubleLF(file_, str);
    element->release();
    serializer->release();
  }
}

void PrsmXmlWriter::writeVector(const PrsmPtrVec &prsm_ptrs) {
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    write(prsm_ptrs[i]);
  }
}

void PrsmXmlWriter::writeVector2D(const PrsmPtrVec2D &prsm_ptrs){
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    writeVector(prsm_ptrs[i]);
  }
}

void PrsmXmlWriter::writeVector3D(const PrsmPtrVec3D &prsm_ptrs){
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    writeVector2D(prsm_ptrs[i]);
  }
}

} /* namespace prot */
