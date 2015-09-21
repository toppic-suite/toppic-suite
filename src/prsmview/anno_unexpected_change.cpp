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

void AnnoUnexpectedChange::addOccurence(int pos,
                                        const std::string &acid_letter) {
  std::pair<int, std::string> new_occurence(pos, acid_letter);
  occurences_.push_back(new_occurence);
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

  if (type_ == "SHIFT") {
      xml_doc->addElement(element, "change_type", convertToString(4).c_str());
  }
  std::string occu = "";
  if (ptm_ptr_ != nullptr) {
      ptm_ptr_->appendxml(xml_doc, element);
      for (size_t i = 0; i < occurences_.size(); i++) {
          if (i != occurences_.size() - 1) {
              occu += occurences_[i].second + convertToString(occurences_[i].first);
              occu += " / ";
          } else {
              occu += occurences_[i].second + convertToString(occurences_[i].first);
          }
      }
  } else {
      if (occurences_.size() > 0) {
          occu += occurences_[0].second + convertToString(occurences_[0].first);
          occu += " - ";
          occu += occurences_[1].second + convertToString(occurences_[1].first);
      }
  }
  xml_doc->addElement(element, "occurence", occu.c_str());
  parent->appendChild(element);
}

}
