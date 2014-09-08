#include "base/anno_change.hpp"

namespace prot {

AnnoChange::AnnoChange(int left_bp_pos, int right_bp_pos, 
                       double mass_shift, int color, int type) {
  left_bp_pos_ = left_bp_pos;
  right_bp_pos_ = right_bp_pos;
  mass_shift_ = mass_shift;
  color_ = color;
  type_ = type;
}

void AnnoChange::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("change");
  std::string str = convertToString(left_bp_pos_);
  xml_doc->addElement(element, "left_bp_pos", str.c_str());
  str = convertToString(right_bp_pos_);
  xml_doc->addElement(element, "right_bp_pos", str.c_str());
  str = convertToString(mass_shift_);
  xml_doc->addElement(element, "mass_shift", str.c_str());
  str = convertToString(color_);
  xml_doc->addElement(element, "unexpected_change_color", str.c_str());
  str = convertToString(type_);
  xml_doc->addElement(element, "type", str.c_str());
  parent->appendChild(element);
}

}
