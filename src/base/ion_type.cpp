#include "base/mass_constant.hpp"
#include "base/ion_type.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

IonType::IonType(const std::string &name, bool n_term, double shift): 
    name_(name),
    n_term_(n_term),
    shift_(shift) {
      if (n_term_) {
        b_y_shift_ = shift_;
      }
      else {
        b_y_shift_ = shift_ - MassConstant::getYIonShift();
      }
    }

IonType::IonType(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  n_term_ = XmlDomUtil::getBoolChildValue(element, "n_term", 0);
  shift_ = XmlDomUtil::getDoubleChildValue(element, "shift", 0);
  if (n_term_) {
    b_y_shift_ = shift_;
  }
  else {
    b_y_shift_ = shift_ - MassConstant::getYIonShift();
  }
}

void IonType::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("ion_type");
  xml_doc->addElement(element, "name", name_.c_str());
  std::string str = StringUtil::convertToString(n_term_);
  xml_doc->addElement(element, "n_term", str.c_str());
  str = StringUtil::convertToString(shift_);
  xml_doc->addElement(element, "shift", str.c_str());
  str = StringUtil::convertToString(b_y_shift_);
  xml_doc->addElement(element, "b_y_shift_", str.c_str());
  parent->appendChild(element);
}

}
