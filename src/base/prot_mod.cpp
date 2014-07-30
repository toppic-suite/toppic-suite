#include <base/logger.hpp>

#include "base/prot_mod.hpp"

namespace prot {

ProtModPtrVec ProtModFactory::prot_mod_ptr_vec_;

ProtMod::ProtMod(const std::string &name, TruncPtr trunc_ptr, 
                 PtmPtr ptm_ptr, const AcidPtrVec &valid_acid_ptr_vec) {
  name_ = name;
  trunc_ptr_ = trunc_ptr;
  ptm_ptr_ = ptm_ptr;
  valid_acid_ptr_vec_ = valid_acid_ptr_vec;
  prot_shift_ = trunc_ptr_->getShift() + ptm_ptr_->getMonoMass();
  pep_shift_ = ptm_ptr_->getMonoMass();
}

void ProtMod::appendxml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("prot_mod");
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

bool ProtMod::allowMod(const ResiduePtrVec &residues){
    if (name_ == "NONE") {
        return true;
    }
    else if (name_ == "NME") {
        if(residues.size()>=2 && residues[0]->getAcidPtr()->getOneLetter() == "M"){
            return true;
        }
        return false;
    }
    else if (name_ == "ACETYLATION") {
        return true;
    }
    else if (name_ == "NME_ACETYLATION"){
        if(residues.size()>=2 && residues[0]->getAcidPtr()->getOneLetter() == "M"){
            return true;
        }
        return false;
    }
  return false;
}

void ProtModFactory::initFactory(const std::string &file_name) {
  prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int mod_num = getChildCount(parent, "prot_mod");
    for (int i = 0; i < mod_num; i++) {
      xercesc::DOMElement* element = getChildElement(parent, "prot_mod", i);
      std::string name = getChildValue(element, "name", 0);
      std::string trunc_name = getChildValue(element, "trunc_name", 0);
      std::string ptm_name = getChildValue(element, "ptm_name", 0);
      std::string valid_acids = getChildValue(element, "valid_acids", 0);
      TruncPtr trunc_ptr = TruncFactory::getBaseTruncPtrByName(trunc_name);
      PtmPtr ptm_ptr = PtmFactory::getBasePtmPtrByAbbrName(ptm_name);
      LOG_DEBUG( "name " << name << " trunc_name " 
                << trunc_name << " valid acids " << valid_acids);
      AcidPtrVec valid_acid_ptrs;
      for (size_t j = 0; j < valid_acids.length(); j++) {
        std::string letter = valid_acids.substr(j, 1);
        AcidPtr acid_ptr = AcidFactory::getBaseAcidPtrByOneLetter(letter);
        valid_acid_ptrs.push_back(acid_ptr);
      }
      prot_mod_ptr_vec_.push_back(ProtModPtr(
              new ProtMod(name, trunc_ptr, ptm_ptr, valid_acid_ptrs)));

    }
  }
}

ProtModPtr ProtModFactory::getBaseProtModPtrByName(const std::string &name) {
  for (size_t i = 0; i < prot_mod_ptr_vec_.size(); i++) {
    std::string n = prot_mod_ptr_vec_[i]->getName();
    if (n == name) {
      return prot_mod_ptr_vec_[i];
    }
  }
  return ProtModPtr(nullptr);
}

}
