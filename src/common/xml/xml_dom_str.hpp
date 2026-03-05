//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_COMMON_XML_XML_DOM_STR_HPP_
#define TOPPIC_COMMON_XML_XML_DOM_STR_HPP_

#include <cassert>
#include <string>

#include <xercesc/util/XMLString.hpp>

namespace toppic {

namespace xml_dom_str {

class XmlStr {
 public:
  explicit XmlStr(const char* str) {
    // Call the private transcoding method
    assert(str != nullptr);
    unicode_form_ = xercesc::XMLString::transcode(str);
  }

  ~XmlStr() {
    if (unicode_form_ != nullptr)
      xercesc::XMLString::release(&unicode_form_);
  }

  const XMLCh* unicodeForm() const {return unicode_form_;}

 private:
  XMLCh* unicode_form_;
};

class CharStr {
 public:
  explicit CharStr(const XMLCh* xml_ch) {
    // Call the private transcoding method
    assert(xml_ch != nullptr);
    ch_ = xercesc::XMLString::transcode(xml_ch);
  }

  ~CharStr() {
    xercesc::XMLString::release(&ch_);
  }

  std::string getString() const {return std::string(ch_);}

 private:
  char* ch_;
};

std::string charStr(const XMLCh* xml_ch);
const XMLCh* xmlStr(const char* str);

}  // namespace xml_dom_str

}  // namespace toppic

#endif  // TOPPIC_COMMON_XML_XML_DOM_STR_HPP_
