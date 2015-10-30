/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <base/logger.hpp>

#include "base/acid.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

Acid::Acid (const std::string &name, const std::string &one_letter, 
            const std::string &three_letter, const std::string &composition, 
            double mono_mass, double avg_mass) {
  name_ = name;
  one_letter_ = one_letter;
  three_letter_ = three_letter;
  composition_ = composition;
  mono_mass_ = mono_mass;
  average_mass_ = avg_mass;
}

Acid::Acid (xercesc::DOMElement* element) { 
  name_ = getChildValue(element, "name", 0);
  one_letter_ = getChildValue(element, "one_letter", 0);
  three_letter_ = getChildValue(element, "three_letter", 0);
  composition_ = getChildValue(element, "composition", 0);
  mono_mass_ = getDoubleChildValue(element, "mono_mass", 0);
  average_mass_ = getDoubleChildValue(element, "average_mass", 0);
}

void Acid::appendxml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = Acid::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  xml_doc->addElement(element, "one_letter", one_letter_.c_str());
  xml_doc->addElement(element, "three_letter", three_letter_.c_str());
  xml_doc->addElement(element, "composition", composition_.c_str());
  std::string str = convertToString(mono_mass_);
  xml_doc->addElement(element, "mono_mass", str.c_str());
  str = convertToString(average_mass_);
  xml_doc->addElement(element, "average_mass", str.c_str());
  parent->appendChild(element);
}

} /* end namespace */

