#include "log4cxx/logger.h"

#include "trunc.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("Trunc"));

Trunc::Trunc(std::string name, int trunc_len, 
                   AcidPtrVec &acid_list, std::string str) {
  name_ = name;
  trunc_len_ = trunc_len;
  shift_ = 0;
  for (unsigned int i = 0; i < str.length(); i++) {
    std::string letter = str.substr(i, 1);
    AcidPtr acid_ptr = getAcidPtrByOneLetter(acid_list, letter);
    acid_str_.push_back(acid_ptr);
    shift_ = shift_ - acid_ptr->getMonoMass();
  }
}


TruncPtrVec getTruncPtrVecInstance(AcidPtrVec &acid_list, 
                                   const char* file_name) {
  TruncPtrVec trunc_list;
  prot::XmlDOMParser* parser = prot::getXmlDOMInstance();
  if (parser) {
    prot::XmlDOMDocument* doc = new prot::XmlDOMDocument(parser, file_name);
    if (doc) {
      xercesc::DOMElement* parent = doc->getDocumentElement();
      int trunc_num = getChildCount(parent, "truncation");
      for (int i = 0; i < trunc_num; i++) {
        xercesc::DOMElement* element = getChildElement(parent, "truncation", i);
        std::string name = getChildValue(element, "name", 0);
        int trunc_len = getIntChildValue(element, "trunc_len", 0);
        std::string str = getChildValue(element, "acid_str", 0);
        LOG4CXX_DEBUG(logger, "name " << name << " str " << str << " trunc len " << trunc_len);
        trunc_list.push_back(TruncPtr(
                new Trunc(name, trunc_len, acid_list, str)));

      }
      delete doc;
    }
  }
  return trunc_list;
}

TruncPtr getTruncPtrByName(TruncPtrVec &trunc_list, 
                         const std::string &name) {
  for (unsigned int i = 0; i < trunc_list.size(); i++) {
    std::string n = trunc_list[i]->getName();
    if (n == name) {
      return trunc_list[i];
    }
  }
  return TruncPtr(nullptr);
}

}
