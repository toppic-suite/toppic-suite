// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
