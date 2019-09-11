//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include <string>

#include <xercesc/util/XMLString.hpp>

namespace toppic {

class XStr {
 public:
  explicit XStr(const char*  str) {
    // Call the private transcoding method
    unicode_form_ = xercesc::XMLString::transcode(str);
  }

  ~XStr() {
    xercesc::XMLString::release(&unicode_form_);
  }

  const XMLCh* unicodeForm() {return unicode_form_;}

 private:
  XMLCh* unicode_form_;
};

class YStr {
 public:
  explicit YStr(const XMLCh* xml_ch) {
    // Call the private transcoding method
    ch_ = xercesc::XMLString::transcode(xml_ch);
  }

  ~YStr() {
    delete ch_;
  }

  std::string  getString() {return std::string(ch_);}

 private:
  char* ch_;
};

#define X(str) XStr(str).unicodeForm()
#define Y(str) YStr(str).getString()

}  // namespace toppic
#endif
