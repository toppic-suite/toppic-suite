// author  Xiaowen Liu
// date    2013-11-1

#include "base/logger.hpp"
#include "base/acid.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Acid::Acid(const std::string &name, const std::string &one_letter, 
           const std::string &three_letter, const std::string &composition, 
           double mono_mass, double average_mass): 
    name_(name),
    one_letter_(one_letter),
    three_letter_(three_letter),
    composition_(composition),
    mono_mass_(mono_mass),
    average_mass_(average_mass) {
    }

Acid::Acid(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  one_letter_ = XmlDomUtil::getChildValue(element, "one_letter", 0);
  three_letter_ = XmlDomUtil::getChildValue(element, "three_letter", 0);
  composition_ = XmlDomUtil::getChildValue(element, "composition", 0);
  mono_mass_ = XmlDomUtil::getDoubleChildValue(element, "mono_mass", 0);
  average_mass_ = XmlDomUtil::getDoubleChildValue(element, "average_mass", 0);
}

void Acid::appendNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = Acid::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

std::string Acid::getNameFromXml(xercesc::DOMElement * element) {
  std::string name = XmlDomUtil::getChildValue(element, "name", 0);
  return name;
}

} 
