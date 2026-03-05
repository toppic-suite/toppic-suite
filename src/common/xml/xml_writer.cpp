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

#include <stdexcept>

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_impl.hpp"
#include "common/xml/xml_dom_util.hpp"

namespace toppic {

XmlWriter::XmlWriter(const std::string &file_name, const std::string &root) {
  file_.open(file_name);
  if (!file_.is_open())
    throw std::runtime_error("XmlWriter: failed to open file: " + file_name);
  root_ = root;
  LOG_DEBUG("file_name " << file_name);
  file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  if (!root_.empty()) {
    file_ << "<" + root_ + ">";
  }
  XmlDOMImpl* impl = XmlDOMImplFactory::getXmlDOMImplInstance();
  doc_ = new XmlDOMDocument(impl->createDoc(!root_.empty() ? root_ : "ROOT"));
  serializer_ = impl->createSerializer();
}

XmlWriter::~XmlWriter() {
  if (file_.is_open()) {
    file_.close();
  }
  serializer_->release();
  delete doc_;
}

// Takes ownership of element and releases it after writing.
void XmlWriter::writeAndRelease(xercesc::DOMElement* element) {
  std::string str = xml_dom_util::writeToString(serializer_, element);
  xml_dom_util::writeToStreamByRemovingDoubleLF(file_, str);
  element->release();
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
