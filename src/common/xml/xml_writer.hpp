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

#include <fstream>
#include <memory>
#include <string>

#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_element.hpp"

namespace toppic {

// Streaming XML writer: it writes the declaration and the opening <root> tag as
// text, then serializes top-level elements one at a time (removing each from the
// scratch document afterwards), so the whole result is never held in memory.
class XmlWriter {
 public:
  XmlWriter(const std::string &file_name, const std::string &root);
  XmlWriter(const XmlWriter&) = delete;
  XmlWriter& operator=(const XmlWriter&) = delete;
  ~XmlWriter();

  // Scratch document used as an element factory. Create a top-level element
  // under getDocumentElement(), build it, then pass it to writeAndRelease().
  XmlDOMDocument* getDoc() { return doc_.get(); }

  // Serialize `element`, write it to the file, then remove it from the scratch
  // document. The element handle must not be used afterwards.
  void writeAndRelease(XmlDOMElement element);

  void writeStr(const std::string& str);

  void close();

 private:
  std::ofstream file_;
  std::string root_;
  // Scratch document: only its root element is used, as the parent under which
  // top-level elements are created before being serialized and removed. The
  // file's root tags are written separately as text, so this is never
  // serialized as a whole.
  std::unique_ptr<XmlDOMDocument> doc_;
};

using XmlWriterPtr = std::shared_ptr<XmlWriter>;

}  // namespace toppic

#endif  // TOPPIC_COMMON_XML_XML_WRITER_HPP_
