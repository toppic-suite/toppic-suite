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

#include <string>

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_parser.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/ion_type_data.hpp"
#include "common/base/ion_type_base.hpp"

namespace toppic {

IonTypePtrVec IonTypeBase::ion_type_ptr_vec_;
IonTypePtr IonTypeBase::ion_type_ptr_B_;
IonTypePtr IonTypeBase::ion_type_ptr_PREC_;

// class functions
void IonTypeBase::initBase() {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing ion type data!");
    exit(EXIT_FAILURE);
  }

  xercesc::MemBufInputSource mem_str((const XMLByte*)ion_type_base_data.c_str(), 
                                     ion_type_base_data.length(), 
                                     "ion_type_data");
  toppic::XmlDOMDocument doc(parser, mem_str);
  XmlDOMElement* parent = doc.getDocumentElement();
  std::string element_name = IonType::getXmlElementName();
  int ion_type_num = xml_dom_util::getChildCount(parent, element_name.c_str());
  for (int i = 0; i < ion_type_num; i++) {
    XmlDOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
    IonTypePtr ion_type_ptr = std::make_shared<IonType>(element);
    ion_type_ptr_vec_.push_back(ion_type_ptr);
    if (ion_type_ptr->getName() == getName_B()) {
      ion_type_ptr_B_ = ion_type_ptr;
    }
    if (ion_type_ptr->getName() == getName_PREC()) {
      ion_type_ptr_PREC_ = ion_type_ptr;
    }
  }
}

IonTypePtr IonTypeBase::getIonTypePtrByName(const std::string &name) {
  for (size_t i = 0; i < ion_type_ptr_vec_.size(); i++) {
    std::string n = ion_type_ptr_vec_[i]->getName();
    if (n == name) {
      return ion_type_ptr_vec_[i];
    }
  }
  LOG_WARN("Ion type " << name << " cannot be found!");
  return IonTypePtr(nullptr);
}

}  // namespace toppic
