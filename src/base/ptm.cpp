/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

#include "base/acid_util.hpp"
#include "base/ptm.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

Ptm::Ptm(const std::string &name, const std::string &abbr_name,
         double mono_mass, int unimod_id, 
         const std::string &n_term_acid_str,
         const std::string &c_term_acid_str, 
         const std::string &anywhere_acid_str): 
    name_(name),
    abbr_name_(abbr_name),
    mono_mass_(mono_mass),
    unimod_id_(unimod_id) {
      n_term_acids_ = AcidUtil::convertStrToAcidPtrVec(n_term_acid_str);
      c_term_acids_ = AcidUtil::convertStrToAcidPtrVec(c_term_acid_str);
      anywhere_acids_ = AcidUtil::convertStrToAcidPtrVec(anywhere_acid_str);
    }

void Ptm::appendNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent) {
  std::string element_name = Ptm::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", abbr_name_.c_str());
  parent->appendChild(element);
}

}

