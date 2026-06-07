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

#include "common/xml/xml_writer.hpp"

#include <memory>
#include <stdexcept>
#include <string>

#include <pugixml.hpp>

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_util.hpp"

namespace toppic {

XmlWriter::XmlWriter(const std::string &file_name, const std::string &root)
    : root_(root) {
  file_.open(file_name);
  if (!file_.is_open()) {
    throw std::runtime_error("XmlWriter: failed to open file: " + file_name);
  }
  LOG_DEBUG("file_name " << file_name);
  file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  if (!root_.empty()) {
    file_ << "<" + root_ + ">";
  }
  // Scratch document whose root is the parent for the top-level elements.
  doc_ = std::make_unique<XmlDOMDocument>(!root_.empty() ? root_ : "ROOT");
}

XmlWriter::~XmlWriter() {
  if (file_.is_open()) {
    file_.close();
  }
  // doc_ is released automatically; pugixml has no serializer to release.
}

void XmlWriter::writeAndRelease(XmlDOMElement element) {
  std::string str = xml_dom_util::writeToString(element);
  xml_dom_util::writeToStreamByRemovingDoubleLF(file_, str);
  // Free the element from the scratch document so the stream stays bounded.
  element.parent().remove_child(element);
}

void XmlWriter::writeStr(const std::string& str) {
  file_ << str << std::endl;
}

void XmlWriter::close() {
  if (!root_.empty()) {
    file_ << "</" + root_ + ">" << std::endl;
  }
  file_.close();
}

}  // namespace toppic
