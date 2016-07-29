#include "base/string_util.hpp"
#include "spec/peak.hpp"

namespace prot {

Peak::Peak(double position, double intensity):
    position_(position),
    intensity_(intensity) {
    }

void Peak::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = Peak::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(getPosition());
  xml_doc->addElement(element, "position", str.c_str());
  str = StringUtil::convertToString(getIntensity());
  xml_doc->addElement(element, "intensity", str.c_str());
  parent->appendChild(element);
}

}
