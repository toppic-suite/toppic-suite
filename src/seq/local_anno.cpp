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
#include <algorithm>

#include "util/str_util.hpp"
#include "xml/xml_dom_util.hpp"
#include "base/ptm_base.hpp"
#include "seq/local_anno.hpp"

namespace toppic {

LocalAnno:: LocalAnno(int left_pos, int right_pos, double conf,
                      const std::vector<double> & scr_vec,
                      double raw_scr,
                      PtmPtr p):
    left_pos_(left_pos),
    right_pos_(right_pos),
    conf_(conf),
    scr_vec_(scr_vec),
    raw_scr_(raw_scr),
    ptm_ptr_(p) {}

LocalAnno::LocalAnno(XmlDOMElement* element) {
  conf_ = xml_dom_util::getDoubleChildValue(element, "confidence", 0);
  std::string scr_str = xml_dom_util::getChildValue(element, "score_list", 0);

  std::vector<std::string> tmp = str_util::split(scr_str, " ");
  for (size_t i = 0; i < tmp.size(); i++) {
    scr_vec_.push_back(std::stod(tmp[i]));
  }

  std::string ptm_element_name = Ptm::getXmlElementName();
  int ptm_count = xml_dom_util::getChildCount(element, ptm_element_name.c_str());

  if (ptm_count == 0) {
    ptm_ptr_ = nullptr;
  } else {
    XmlDOMElement* ptm_element
        = xml_dom_util::getChildElement(element, ptm_element_name.c_str(), 0);
    ptm_ptr_ = PtmBase::getPtmPtrFromXml(ptm_element);
  }
}

double LocalAnno::getScr() {
  return std::accumulate(scr_vec_.begin(), scr_vec_.end(), 0.0);
}

void LocalAnno::appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = str_util::toString(conf_, 4);
  xml_doc->addElement(element, "confidence", str.c_str());

  str = str_util::toString(scr_vec_[0], 4);
  for (size_t i = 1; i < scr_vec_.size(); i++) {
    str = str + " " + str_util::toString(scr_vec_[i], 4);
  }

  xml_doc->addElement(element, "score_list", str.c_str());

  if (ptm_ptr_ != nullptr) {
    ptm_ptr_->appendAbbrNameToXml(xml_doc, element);
  }
  parent->appendChild(element);
}

}  // namespace toppic
