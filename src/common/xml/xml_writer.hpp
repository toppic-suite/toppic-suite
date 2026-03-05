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

#ifndef TOPPIC_COMMON_XML_XML_WRITER_HPP_
#define TOPPIC_COMMON_XML_XML_WRITER_HPP_

#include <string>
#include <fstream>
#include <memory>

#include "common/xml/xml_dom_document.hpp"

namespace toppic {

class XmlWriter {
 public:
  XmlWriter(const std::string &file_name,
            const std::string &root);
  XmlWriter(const XmlWriter&) = delete;
  XmlWriter& operator=(const XmlWriter&) = delete;

  ~XmlWriter();

  XmlDOMDocument* getDoc() { return doc_; }

  void writeAndRelease(xercesc::DOMElement* element);

  void writeStr(const std::string& str);

  void close();

 private:
  xercesc::DOMLSSerializer* serializer_;

  XmlDOMDocument* doc_;

  std::ofstream file_;

  std::string root_;
};

using XmlWriterPtr = std::shared_ptr<XmlWriter>;

}  // namespace toppic

#endif  // TOPPIC_COMMON_XML_XML_WRITER_HPP_
