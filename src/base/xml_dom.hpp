#ifndef PROT_XML_DOM_HPP_
#define PROT_XML_DOM_HPP_

#include <string>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>

namespace prot {

/* DOM parser */
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

 
class XStr {
 public:
  XStr(const char*  str) {
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
  YStr(const XMLCh* xml_ch) {
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

}
#endif
