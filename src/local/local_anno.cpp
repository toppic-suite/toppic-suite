// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "base/ptm_base.hpp"
#include "local_anno.hpp"

namespace prot {

LocalAnno::LocalAnno(xercesc::DOMElement* element) {
  conf_ = XmlDomUtil::getDoubleChildValue(element, "confidence", 0);
  std::string scr_str = XmlDomUtil::getChildValue(element, "score_list", 0);
  std::vector<std::string> tmp = StringUtil::split(scr_str, ' ');
  for (size_t i = 0; i < tmp.size(); i++) {
    scr_vec_.push_back(std::stod(tmp[i]));
  }
  std::string ptm_element_name = Ptm::getXmlElementName();
  int ptm_count = XmlDomUtil::getChildCount(element, ptm_element_name.c_str());
  if (ptm_count == 0) {
    ptm_ptr_ = nullptr;
  } else {
    xercesc::DOMElement* ptm_element 
        = XmlDomUtil::getChildElement(element, ptm_element_name.c_str(), 0);
    ptm_ptr_ = PtmBase::getPtmPtrFromXml(ptm_element);        
  }
}

void LocalAnno::appendToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent) {
  std::string element_name = getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(conf_);
  xml_doc->addElement(element, "confidence", str.c_str());

  str = StringUtil::convertToString(scr_vec_[0]);
  for (size_t i = 1; i < scr_vec_.size(); i++) {
    str = str + " " + StringUtil::convertToString(scr_vec_[i]);
  }

  xml_doc->addElement(element, "score_list", str.c_str());

  if (ptm_ptr_ != nullptr) {
    ptm_ptr_->appendAbbrNameToXml(xml_doc, element);
  }
  parent->appendChild(element);
}

}
