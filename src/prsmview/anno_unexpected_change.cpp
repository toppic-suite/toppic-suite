#include "prsmview/anno_unexpected_change.hpp"

namespace prot {

AnnoUnexpectedChange::AnnoUnexpectedChange(int left_pos, int right_pos, 
                                           double mass_shift, int color, 
                                           const std::string &type) {
  left_pos_ = left_pos;
  right_pos_ = right_pos;
  mass_shift_ = mass_shift;
  color_ = color;
  type_ = type;
}

void AnnoUnexpectedChange::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent,
                                     int decimal_point_num){
  xercesc::DOMElement* element = xml_doc->createElement("unexpected_change");
  std::string str = convertToString(left_pos_);
  xml_doc->addElement(element, "left_position", str.c_str());
  str = convertToString(right_pos_);
  xml_doc->addElement(element, "right_position", str.c_str());
  str = convertToString(mass_shift_, decimal_point_num);
  xml_doc->addElement(element, "mass_shift", str.c_str());
  str = convertToString(color_);
  xml_doc->addElement(element, "unexpected_change_color", str.c_str());
  xml_doc->addElement(element, "type", type_.c_str());
  parent->appendChild(element);
}

}
