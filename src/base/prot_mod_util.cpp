#include "base/logger.hpp"
#include "base/mod_base.hpp"
#include "base/prot_mod_base.hpp"
#include "base/prot_mod_util.hpp"
#include "base/trunc_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

bool ProtModUtil::allowMod(ProtModPtr prot_mod_ptr, const ResiduePtrVec &residues){
  if (prot_mod_ptr == ProtModBase::getProtModPtr_NONE()) {
    return true;
  }
  else {
    // check trunc
    if (!TruncUtil::isValidTrunc(prot_mod_ptr->getTruncPtr(), residues)) {
      return false;
    }
    ModPtr mod_ptr = prot_mod_ptr->getModPtr();
    if (mod_ptr != ModBase::getNoneModPtr()) {
      int mod_pos = prot_mod_ptr->getModPos();
      if (mod_pos < (int)residues.size()) { 
        return false;
      }
      if (residues[mod_pos] != mod_ptr->getOriResiduePtr()) {
        return false;
      }
    }
    return true;
  }
}

ProtModPtrVec ProtModUtil::readProtMod(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  ProtModPtrVec mod_ptr_vec;
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = ProtMod::getXmlElementName();
    int mod_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    LOG_DEBUG("mod num " << mod_num);
    for (int i = 0; i < mod_num; i++) {
      xercesc::DOMElement* element 
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ProtModPtr ptr = ProtModBase::getProtModPtrFromXml(element);
      mod_ptr_vec.push_back(ptr);
    }
  }
  return mod_ptr_vec;
}


/*
bool ProtModUtil::contain_NME_ACETYLATION(const ProtModPtrVec &prot_mod_ptrs) {
  for (size_t i = 0; i < prot_mod_ptrs.size(); i++) {
    if (prot_mod_ptrs[i] == ProtModBase::getProtModPtr_NME_ACETYLATION()) {
      return true;
    }
  }
  return false;
}
*/

}
