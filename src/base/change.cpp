#include "base/ptm_base.hpp"
#include "base/change.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"

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
  left_bp_pos_ = XmlDomUtil::getIntChildValue(element, "left_bp_pos", 0);
  right_bp_pos_ = XmlDomUtil::getIntChildValue(element, "right_bp_pos", 0);
  std::string ct_element_name = ChangeType::getXmlElementName();
  xercesc::DOMElement* ct_element 
      = XmlDomUtil::getChildElement(element, ct_element_name.c_str(), 0);
  change_type_ptr_ = ChangeType::getChangeTypePtrFromXml(ct_element);
  mass_shift_ = XmlDomUtil::getDoubleChildValue(element, "mass_shift", 0);
  std::string ptm_element_name = Ptm::getXmlElementName();
  int ptm_count = XmlDomUtil::getChildCount(element, ptm_element_name.c_str());
  if (ptm_count != 0) {
    xercesc::DOMElement* ptm_element 
        = XmlDomUtil::getChildElement(element, ptm_element_name.c_str(), 0);
    ptm_ptr_ = PtmBase::getPtmPtrFromXml(ptm_element);
  }
}

void Change::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = Change::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(left_bp_pos_);
  xml_doc->addElement(element, "left_bp_pos", str.c_str());
  str = StringUtil::convertToString(right_bp_pos_);
  xml_doc->addElement(element, "right_bp_pos", str.c_str());
  change_type_ptr_->appendXml(xml_doc, element);
  str = StringUtil::convertToString(mass_shift_);
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

ChangePtr Change::geneChangePtr(ChangePtr ori_ptr, int start_pos) {
  int left_bp_pos = ori_ptr->left_bp_pos_ - start_pos;
  int right_bp_pos = ori_ptr->right_bp_pos_ - start_pos;
  ChangeTypePtr change_type_ptr = ori_ptr->change_type_ptr_;
  double mass_shift = ori_ptr->mass_shift_;
  PtmPtr ptm_ptr = ori_ptr->ptm_ptr_;
  ChangePtr change_ptr(
      new Change(left_bp_pos, right_bp_pos, change_type_ptr, mass_shift, ptm_ptr));
  return change_ptr;
}

}
