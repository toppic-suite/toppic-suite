#include <iostream>

#include "acid_util.hpp"
#include "ptm_util.hpp"
#include "trunc_util.hpp"
#include "prot_mod_util.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace prot {

ProtModPtrVec getProtModPtrVecInstance(AcidPtrVec &acid_ptr_vec,
                                       PtmPtrVec &ptm_ptr_vec,
                                       TruncPtrVec &trunc_ptr_vec,
                                       const char* file_name) {
  ProtModPtrVec prot_mod_ptr_vec;
  prot::XmlDOMParser* parser = prot::getXmlDOMInstance();
  if (parser) {
    prot::XmlDOMDocument* doc = new prot::XmlDOMDocument(parser, file_name);
    if (doc) {
      int mod_num = doc->getChildCount("mod_list", 0, "mod");
      for (int i = 0; i < mod_num; i++) {
        xercesc::DOMElement* element = doc->getElement("mod", i);
        std::string name = getChildValue(element, "name");
        std::string trunc_name = getChildValue(element, "trunc_name");
        std::string ptm_name = getChildValue(element, "ptm_name");
        std::string valid_acids = getChildValue(element, "valid_acid");
        TruncPtr trunc_ptr = getTruncPtrByName(trunc_ptr_vec, trunc_name);
        PtmPtr ptm_ptr = getPtmPtrByAbbrName(ptm_ptr_vec, ptm_name);
        AcidPtrVec valid_acid_ptrs;
        for (unsigned int j = 0; j < valid_acids.length(); j++) {
          std::string letter = valid_acids.substr(j, 1);
          AcidPtr acid_ptr = getAcidPtrByOneLetter(acid_ptr_vec, letter);
          valid_acid_ptrs.push_back(acid_ptr);
        }
        prot_mod_ptr_vec.push_back(ProtModPtr(
                new ProtMod(name, trunc_ptr, ptm_ptr, valid_acid_ptrs)));

      }
      delete doc;
    }
    delete parser;
  }
  return prot_mod_ptr_vec;
}

ProtModPtr getProtModPtrByName(ProtModPtrVec &prot_mod_ptr_vec, 
                         const std::string &name) {
  for (unsigned int i = 0; i < prot_mod_ptr_vec.size(); i++) {
    std::string n = prot_mod_ptr_vec[i]->getName();
    if (n.compare(name) == 0) {
      return prot_mod_ptr_vec[i];
    }
  }
  return ProtModPtr(nullptr);
}

}
