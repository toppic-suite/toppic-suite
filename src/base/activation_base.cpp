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


#include "base/activation_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ActivationPtrVec ActivationBase::activation_ptr_vec_;

void ActivationBase::initBase(const std::string &file_name){
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Activation::getXmlElementName();
    int activation_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < activation_num; i++) {
      xercesc::DOMElement* element 
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ActivationPtr ptr(new Activation(element));
      activation_ptr_vec_.push_back(ptr);
    }
  }
}

ActivationPtr ActivationBase::getActivationPtrByName(
    const std::string &name){
  for (size_t i = 0; i < activation_ptr_vec_.size(); i++) {
    std::string n = activation_ptr_vec_[i]->getName();
    if (n == name) {
      return activation_ptr_vec_[i];
    }
  }
  return ActivationPtr(nullptr);
}

ActivationPtr ActivationBase::getActivationPtrFromXml(xercesc::DOMElement * element) {
  std::string name = Activation::getNameFromXml(element);
  ActivationPtr activation_ptr = getActivationPtrByName(name);
  return activation_ptr;
}

} 
