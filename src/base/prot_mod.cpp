#include "prot_mod.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace prot {

ProtMod::ProtMod(std::string name, TruncPtr trunc_ptr, PtmPtr ptm_ptr,
                 AcidPtrVec valid_acid_ptrs) {
  name_ = name;
  trunc_ptr_ = trunc_ptr;
  ptm_ptr_ = ptm_ptr;
  valid_acid_ptrs_ = valid_acid_ptrs;
  prot_shift_ = trunc_ptr_->getShift() + ptm_ptr_->getMonoMass();
  pep_shift_ = ptm_ptr_->getMonoMass();
}

ProtModPtrVec getProtModPtrVecInstance(AcidPtrVec &acid_list,
                                       PtmPtrVec &ptm_list,
                                       TruncPtrVec &trunc_list,
                                       const char* file_name) {
  ProtModPtrVec prot_mod_list;
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
        TruncPtr trunc_ptr = getTruncPtrByName(trunc_list, trunc_name);
        PtmPtr ptm_ptr = getPtmPtrByAbbrName(ptm_list, ptm_name);
        AcidPtrVec valid_acid_ptrs;
        for (unsigned int j = 0; j < valid_acids.length(); j++) {
          std::string letter = valid_acids.substr(j, 1);
          AcidPtr acid_ptr = getAcidPtrByOneLetter(acid_list, letter);
          valid_acid_ptrs.push_back(acid_ptr);
        }
        prot_mod_list.push_back(ProtModPtr(
                new ProtMod(name, trunc_ptr, ptm_ptr, valid_acid_ptrs)));

      }
      delete doc;
    }
    delete parser;
  }
  return prot_mod_list;
}

ProtModPtr getProtModPtrByName(ProtModPtrVec &prot_mod_list, 
                         const std::string &name) {
  for (unsigned int i = 0; i < prot_mod_list.size(); i++) {
    std::string n = prot_mod_list[i]->getName();
    if (n == name) {
      return prot_mod_list[i];
    }
  }
  return ProtModPtr(nullptr);
}

}
