/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

#include "base/ptm.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Ptm::Ptm(const std::string &name, const std::string &abbr_name,
         double mono_mass, int unimod_id, 
         const std::string &n_term_residue_str,
         const std::string &c_term_residue_str, 
         const std::string &anywhere_residue_str): 
    name_(name),
    abbr_name_(abbr_name),
    mono_mass_(mono_mass),
    unimod_id_(unimod_id) {
    }

Ptm::Ptm(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  abbr_name_ = XmlDomUtil::getChildValue(element, "abbreviation", 0);
  mono_mass_ = XmlDomUtil::getDoubleChildValue(element, "mono_mass", 0);
  unimod_id_ = XmlDomUtil::getIntChildValue(element, "unimod", 0);
}

void Ptm::appendAbbrNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent) {
  std::string element_name = Ptm::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "abbreviation", abbr_name_.c_str());
  parent->appendChild(element);
}

std::string Ptm::getAbbrNameFromXml(xercesc::DOMElement * element) {
  std::string abbr_name = XmlDomUtil::getChildValue(element, "abbreviation", 0);
  return abbr_name;
}

}

