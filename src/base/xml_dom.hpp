#ifndef PROT_XML_DOM_HPP_
#define PROT_XML_DOM_HPP_

#include <string>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/sax/HandlerBase.hpp>

namespace prot {

class XmlDOMParser {
 public:
  XmlDOMParser();
  ~XmlDOMParser();

  xercesc::DOMDocument* parse(const char* xml_file);

 private:
  xercesc::XercesDOMParser* parser_;
  xercesc::ErrorHandler*    err_handler_;

};

class XmlDOMParserFactory {
 private:
  static XmlDOMParser* dom_parser;
 public:
  static XmlDOMParser* getXmlDOMInstance() {
    if (dom_parser == nullptr) {
      dom_parser = new XmlDOMParser();
    }
    return dom_parser;
  }
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
    xercesc::XMLString::release(&ch_);
  }

  const char * getString() {return ch_;}

 private:
  char * ch_;
};

#define X(str) XStr(str).unicodeForm()
#define Y(str) YStr(str).getString()

int writeXmlFile(xercesc::DOMDocument* doc, const char * filename);

}
#endif
