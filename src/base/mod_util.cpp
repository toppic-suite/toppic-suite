#include "base/logger.hpp"
#include "base/mod_base.hpp"
#include "base/mod_util.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ModPtrVec ModUtil::readMod(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  ModPtrVec mod_ptr_vec;
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Mod::getXmlElementName();
    int mod_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    LOG_DEBUG("mod num " << mod_num);
    for (int i = 0; i < mod_num; i++) {
      xercesc::DOMElement* element 
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ModPtr ptr = ModBase::getModPtrFromXml(element);
      mod_ptr_vec.push_back(ptr);
    }
  }
  return mod_ptr_vec;
}

ModPtrVec ModUtil::geneFixedModList(const std::string &str) {
  if (str == "" || str == "C57" || str == "C58") {
    ModPtrVec mod_ptr_vec;
    if (str == "C57") {
      mod_ptr_vec.push_back(ModBase::getC57ModPtr());
    }
    else if (str == "C58") {
      mod_ptr_vec.push_back(ModBase::getC58ModPtr());
    }
    return mod_ptr_vec;
  }
  else {
    return readMod(str);
  }
}

ResiduePtrVec ModUtil::geneResidueListWithMod(ResiduePtrVec residue_list,
                                              ModPtrVec fix_mod_list) {
  ResiduePtrVec result;
  for (size_t i = 0; i < residue_list.size(); i++) {
    bool mod = false;
    for (size_t j = 0; j < fix_mod_list.size(); j++) {
      if (fix_mod_list[j]->getOriResiduePtr() == residue_list[i]) {
        mod = true;
        result.push_back(fix_mod_list[j]->getModResiduePtr());
        break;                 
      }
    }
    if (!mod) {
      result.push_back(residue_list[i]);
    }
  }
  return result;
}

} /* end namespace */

