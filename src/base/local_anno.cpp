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


#include <string>
#include <vector>

#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "base/ptm_base.hpp"
#include "base/local_anno.hpp"

namespace toppic {

LocalAnno::LocalAnno(xercesc::DOMElement* element) {
  conf_ = xml_dom_util::getDoubleChildValue(element, "confidence", 0);
  std::string scr_str = xml_dom_util::getChildValue(element, "score_list", 0);

  std::vector<std::string> tmp = string_util::split(scr_str, ' ');
  for (size_t i = 0; i < tmp.size(); i++) {
    scr_vec_.push_back(std::stod(tmp[i]));
  }

  std::string ptm_element_name = Ptm::getXmlElementName();
  int ptm_count = xml_dom_util::getChildCount(element, ptm_element_name.c_str());

  if (ptm_count == 0) {
    ptm_ptr_ = nullptr;
  } else {
    xercesc::DOMElement* ptm_element
        = xml_dom_util::getChildElement(element, ptm_element_name.c_str(), 0);
    ptm_ptr_ = PtmBase::getPtmPtrFromXml(ptm_element);
  }
}

void LocalAnno::appendToXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = string_util::convertToString(conf_, 4);
  xml_doc->addElement(element, "confidence", str.c_str());

  str = string_util::convertToString(scr_vec_[0], 4);
  for (size_t i = 1; i < scr_vec_.size(); i++) {
    str = str + " " + string_util::convertToString(scr_vec_[i], 4);
  }

  xml_doc->addElement(element, "score_list", str.c_str());

  if (ptm_ptr_ != nullptr) {
    ptm_ptr_->appendAbbrNameToXml(xml_doc, element);
  }
  parent->appendChild(element);
}

}  // namespace toppic
