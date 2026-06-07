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

#include <stdexcept>
#include <string>

#include <xercesc/util/XMLString.hpp>

namespace toppic {

class XmlStr {
 public:
  explicit XmlStr(const char* str) {
    if (str == nullptr)
      throw std::invalid_argument("XmlStr: input str is null");
    unicode_form_ = xercesc::XMLString::transcode(str);
    if (unicode_form_ == nullptr)
      throw std::runtime_error("XmlStr: XMLString::transcode failed");
  }

  XmlStr(const XmlStr&) = delete;
  XmlStr& operator=(const XmlStr&) = delete;

  ~XmlStr() {
    if (unicode_form_ != nullptr)
      xercesc::XMLString::release(&unicode_form_);
  }

  const XMLCh* unicodeForm() const { return unicode_form_; }

 private:
  XMLCh* unicode_form_;
};

class CharStr {
 public:
  explicit CharStr(const XMLCh* xml_ch) {
    if (xml_ch == nullptr)
      throw std::invalid_argument("CharStr: input xml_ch is null");
    ch_ = xercesc::XMLString::transcode(xml_ch);
    if (ch_ == nullptr)
      throw std::runtime_error("CharStr: XMLString::transcode failed");
  }

  CharStr(const CharStr&) = delete;
  CharStr& operator=(const CharStr&) = delete;

  ~CharStr() {
    if (ch_ != nullptr)
      xercesc::XMLString::release(&ch_);
  }

  const char* c_str() const { return ch_; }
  std::string getString() const { return std::string(ch_); }

 private:
  char* ch_;
};

}  // namespace toppic

#endif  // TOPPIC_COMMON_XML_XML_DOM_STR_HPP_
