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

#include "common/xml/xml_dom_impl.hpp"

#include <stdexcept>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMLSSerializer.hpp>

#include "common/xml/xml_dom_str.hpp"

namespace toppic {

/* XmlDOMImplementation */
XmlDOMImpl* XmlDOMImplFactory::dom_impl_ = nullptr;

XmlDOMImpl::XmlDOMImpl() {
  impl_ = xercesc::DOMImplementationRegistry::getDOMImplementation(XmlStr("Core").unicodeForm());
  if (impl_ == nullptr)
    throw std::runtime_error("XmlDOMImpl: DOMImplementation unavailable");
}

xercesc::DOMDocument* XmlDOMImpl::createDoc(const std::string &root) {
  return impl_->createDocument(0, XmlStr(root.c_str()).unicodeForm(), 0);
}

xercesc::DOMLSSerializer* XmlDOMImpl::createSerializer() {
  xercesc::DOMLSSerializer* writer = impl_->createLSSerializer();
  writer->getDomConfig()->setParameter(
      xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);
  writer->getDomConfig()->setParameter(
      xercesc::XMLUni::fgDOMWRTDiscardDefaultContent, true);
  writer->setNewLine(XmlStr("\n").unicodeForm());
  return writer;
}

/* XmlDOMImplFactory */
XmlDOMImpl* XmlDOMImplFactory::getXmlDOMImplInstance() {
  if (dom_impl_ == nullptr) {
    dom_impl_ = new XmlDOMImpl();
  }
  return dom_impl_;
}

}  // namespace toppic
