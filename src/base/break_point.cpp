
#include "base/break_point.hpp"

namespace prot {

BreakPoint::BreakPoint(double prm, double srm){
  prm_ = prm;
  srm_ = srm;
}

void BreakPoint::appendXml(XmlDOMDocument* xml_doc,
                           xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("breakpoint");
  std::string str = convertToString(prm_);
  xml_doc->addElement(element, "prm", str.c_str());
  str = convertToString(srm_);
  xml_doc->addElement(element, "srm", str.c_str());
  parent->appendChild(element);
}

} /* namespace prot */
