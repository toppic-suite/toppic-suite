#include <iostream>

#include "trunc_util.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace prot {

TruncPtrVec getTruncPtrVecInstance(AcidPtrVec &acid_ptr_vec, 
                                         const char* file_name) {
  TruncPtrVec trunc_ptr_vec;
  prot::XmlDOMParser* parser = prot::getXmlDOMInstance();
  if (parser) {
    prot::XmlDOMDocument* doc = new prot::XmlDOMDocument(parser, file_name);
    if (doc) {
      int acid_num = doc->getChildCount("exp_trunc_list", 0, "exp_trunc");
      for (int i = 0; i < acid_num; i++) {
        xercesc::DOMElement* element = doc->getElement("amino_acid", i);
        std::string name = getChildValue(element, "name");
        int trunc_len = getIntChildValue(element, "trunc_len");
        std::string acid_str = getChildValue(element, "acid_str");
        trunc_ptr_vec.push_back(TruncPtr(
                new Trunc(name, trunc_len, acid_ptr_vec, acid_str)));

      }
      delete doc;
    }
    delete parser;
  }
  return trunc_ptr_vec;
}

TruncPtr getTruncPtrByName(TruncPtrVec &trunc_ptr_vec, 
                         const std::string &name) {
  for (unsigned int i = 0; i < trunc_ptr_vec.size(); i++) {
    std::string n = trunc_ptr_vec[i]->getName();
    if (n.compare(name) == 0) {
      return trunc_ptr_vec[i];
    }
  }
  return TruncPtr(nullptr);
}

}


