#include "base/xml_dom_util.hpp"
#include "base/string_util.hpp"
#include "prsm/n_term_shift.hpp"

namespace prot {

NTermShift:: NTermShift(double shift, int type):
    shift_(shift),
    type_(type) {
    }

NTermShift::NTermShift(xercesc::DOMElement* element){
  shift_ = XmlDomUtil::getDoubleChildValue(element, "shift", 0);
  type_ = XmlDomUtil::getIntChildValue(element, "type", 0);
}

void NTermShift::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent){
  std::string element_name = NTermShift::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(shift_);
  xml_doc->addElement(element, "shift", str.c_str());
  str = StringUtil::convertToString(type_);
  xml_doc->addElement(element, "type", str.c_str());
  parent->appendChild(element);
}

} /* namespace prot */
