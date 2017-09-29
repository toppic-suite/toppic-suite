// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_BASE_XML_DOM_HPP_
#define PROT_BASE_XML_DOM_HPP_

#include <string>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>

namespace prot {

// DOM parser
class XmlDOMParser {
 public:
  XmlDOMParser();
  ~XmlDOMParser();

  xercesc::DOMDocument* parse(const std::string &xml_file);

  xercesc::DOMDocument* parse(const xercesc::MemBufInputSource &str_buf);

 private:
  xercesc::XercesDOMParser* parser_;
  xercesc::ErrorHandler*    err_handler_;

};

class XmlDOMParserFactory {
 private:
  static XmlDOMParser* dom_parser_;
 public:
  static XmlDOMParser* getXmlDOMParserInstance(); 
};

/* DOM Implementation */
class XmlDOMImpl{
 public:
  XmlDOMImpl();
  ~XmlDOMImpl();
  xercesc::DOMDocument* createDoc(const std::string &root);
  xercesc::DOMLSSerializer* createSerializer();

 private:
  xercesc::DOMImplementation* impl_;
};

class XmlDOMImplFactory {
 private:
  static XmlDOMImpl* dom_impl_;
 public:
  static XmlDOMImpl* getXmlDOMImplInstance();
};

inline const XMLCh* X(const char*  str) {
  return xercesc::XMLString::transcode(str);
}

inline std::string Y(const XMLCh* xml_ch) {
  return std::string(xercesc::XMLString::transcode(xml_ch));
}
}  // namespace prot
#endif
