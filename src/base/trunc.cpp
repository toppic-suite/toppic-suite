#include <base/logger.hpp>

#include "base/trunc.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

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
                                   const std::string &file_name) {
  TruncPtrVec trunc_list;
  prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int trunc_num = getChildCount(parent, "truncation");
    for (int i = 0; i < trunc_num; i++) {
      xercesc::DOMElement* element = getChildElement(parent, "truncation", i);
      std::string name = getChildValue(element, "name", 0);
      int trunc_len = getIntChildValue(element, "trunc_len", 0);
      std::string str = getChildValue(element, "acid_str", 0);
      LOG_DEBUG( "name " << name << " str " << str << " trunc len " << trunc_len);
      trunc_list.push_back(TruncPtr(
              new Trunc(name, trunc_len, acid_list, str)));
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
