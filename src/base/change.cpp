#include "base/ptm_base.hpp"
#include "base/change.hpp"

namespace prot {

Change::Change(int left_bp_pos, int right_bp_pos, 
               ChangeTypePtr change_type_ptr,
               double mass_shift, PtmPtr ptm_ptr): 
    left_bp_pos_(left_bp_pos), 
    right_bp_pos_(right_bp_pos),
    change_type_ptr_(change_type_ptr),
    mass_shift_(mass_shift),
    ptm_ptr_(ptm_ptr) {
    }

Change::Change(xercesc::DOMElement* element) {
  left_bp_pos_ = getIntChildValue(element, "left_bp_pos", 0);
  right_bp_pos_ = getIntChildValue(element, "right_bp_pos", 0);
  std::string ct_element_name = ChangeType::getXmlElementName();
  xercesc::DOMElement* ct_element 
      = getChildElement(element, ct_element_name.c_str(), 0);
  change_type_ptr_ = ChangeType::getChangeTypePtrFromXml(ct_element);
  mass_shift_ = getDoubleChildValue(element, "mass_shift", 0);
  std::string ptm_element_name = Ptm::getXmlElementName();
  int ptm_count = getChildCount(element, ptm_element_name.c_str());
  if (ptm_count != 0) {
    xercesc::DOMElement* ptm_element 
        = getChildElement(element, ptm_element_name.c_str(), 0);
    ptm_ptr_ = PtmBase::getPtmPtrFromXml(ptm_element);
  }
}

void Change::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = Change::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = convertToString(left_bp_pos_);
  xml_doc->addElement(element, "left_bp_pos", str.c_str());
  str = convertToString(right_bp_pos_);
  xml_doc->addElement(element, "right_bp_pos", str.c_str());
  change_type_ptr_->appendXml(xml_doc, element);
  str = convertToString(mass_shift_);
  xml_doc->addElement(element, "mass_shift", str.c_str());
  if (ptm_ptr_ != nullptr) {
    ptm_ptr_->appendAbbrNameToXml(xml_doc, element);
  }
  parent->appendChild(element);
}

bool Change::cmpPosIncrease(const ChangePtr &a, const ChangePtr &b) {
  if (a->getLeftBpPos() < b->getLeftBpPos()) {
    return true;
  }
  else if (a->getLeftBpPos() > b->getLeftBpPos()) {
    return false;
  }
  else {
    return a->getRightBpPos() < b->getRightBpPos();
  }
}

}
