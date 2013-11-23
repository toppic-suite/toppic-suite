#include <iostream>

#include "exp_trunc_util.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace prot {

ExpTruncPtrVec getExpTruncPtrVecInstance(AcidPtrVec &acid_ptr_vec, 
                                         const char* file_name) {
  ExpTruncPtrVec exp_trunc_ptr_vec;
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
        exp_trunc_ptr_vec.push_back(ExpTruncPtr(
                new ExpTrunc(name, trunc_len, acid_ptr_vec, acid_str)));

      }
      delete doc;
    }
    delete parser;
  }
  return exp_trunc_ptr_vec;
}

ExpTruncPtr getExpTruncPtrByName(ExpTruncPtrVec &exp_trunc_ptr_vec, 
                         const std::string &name) {
  for (unsigned int i = 0; i < exp_trunc_ptr_vec.size(); i++) {
    std::string n = exp_trunc_ptr_vec[i]->getName();
    if (n.compare(name) == 0) {
      return exp_trunc_ptr_vec[i];
    }
  }
  return ExpTruncPtr(nullptr);
}

}


