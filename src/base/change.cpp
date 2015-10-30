#include "base/change.hpp"

namespace prot {

Change::Change(int left_bp_pos, int right_bp_pos, int change_type,
               double mass_shift, const PtmPtr & ptm_ptr) {
  left_bp_pos_ = left_bp_pos;
  right_bp_pos_ = right_bp_pos;
  change_type_ = change_type;
  mass_shift_ = mass_shift;
  ptm_ptr_ = ptm_ptr;
}
  left_bp_pos_ = getIntChildValue(change_element, "left_bp_pos", 0);
  right_bp_pos_ = getIntChildValue(change_element, "right_bp_pos", 0);
  change_type_ = getIntChildValue(change_element, "change_type", 0);
  mass_shift_ = getDoubleChildValue(change_element, "mass_shift", 0);
  int ptm_count = getChildCount(change_element, "ptm");
  if (ptm_count != 0) {
    xercesc::DOMElement* ptm_element 
        = getChildElement(change_element, "ptm", 0);
    std::string ptm_abbr_name = getChileElement(ptm_element, "abbr_name", 0);
    ptm_ptr_ = PtmFactory::getBasePtmPtrByAbbrName(ptm_abbr_name);
  }
}

void Change::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("change");
  std::string str = convertToString(left_bp_pos_);
  xml_doc->addElement(element, "left_bp_pos", str.c_str());
  str = convertToString(right_bp_pos_);
  xml_doc->addElement(element, "right_bp_pos", str.c_str());
  str = convertToString(change_type_);
  xml_doc->addElement(element, "change_type", str.c_str());
  str = convertToString(mass_shift_);
  xml_doc->addElement(element, "mass_shift", str.c_str());

  if (ptm_ptr_ != nullptr) {
    ptm_ptr_->appendxml(xml_doc, element);
  }
  parent->appendChild(element);
}

}
