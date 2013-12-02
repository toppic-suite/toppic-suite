#include <log4cxx/logger.h>

#include "base/prot_mod.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("ProtMod"));

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
      xercesc::DOMElement* parent = doc->getDocumentElement();
      int mod_num = getChildCount(parent, "prot_mod");
      for (int i = 0; i < mod_num; i++) {
        xercesc::DOMElement* element = getChildElement(parent, "prot_mod", i);
        std::string name = getChildValue(element, "name", 0);
        std::string trunc_name = getChildValue(element, "trunc_name", 0);
        std::string ptm_name = getChildValue(element, "ptm_name", 0);
        std::string valid_acids = getChildValue(element, "valid_acids", 0);
        TruncPtr trunc_ptr = getTruncPtrByName(trunc_list, trunc_name);
        PtmPtr ptm_ptr = getPtmPtrByAbbrName(ptm_list, ptm_name);
        LOG4CXX_DEBUG(logger, "name " << name << " trunc_name " << trunc_name << " valid acids " << valid_acids);
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
