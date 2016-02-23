#include "base/mod_base.hpp"
#include "base/change.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Change::Change(int left_bp_pos, int right_bp_pos, 
               ChangeTypePtr change_type_ptr,
               double mass_shift, ModPtr mod_ptr): 
    left_bp_pos_(left_bp_pos), 
    right_bp_pos_(right_bp_pos),
    change_type_ptr_(change_type_ptr),
    mass_shift_(mass_shift),
    mod_ptr_(mod_ptr) {
      local_anno_ptr_ = nullptr;
    }

Change::Change(xercesc::DOMElement* element) {
  left_bp_pos_ = XmlDomUtil::getIntChildValue(element, "left_bp_pos", 0);
  right_bp_pos_ = XmlDomUtil::getIntChildValue(element, "right_bp_pos", 0);
  std::string ct_element_name = ChangeType::getXmlElementName();
  xercesc::DOMElement* ct_element 
      = XmlDomUtil::getChildElement(element, ct_element_name.c_str(), 0);
  change_type_ptr_ = ChangeType::getChangeTypePtrFromXml(ct_element);
  mass_shift_ = XmlDomUtil::getDoubleChildValue(element, "mass_shift", 0);
  std::string mod_element_name = Mod::getXmlElementName();
  int mod_count = XmlDomUtil::getChildCount(element, mod_element_name.c_str());
  if (mod_count != 0) {
    xercesc::DOMElement* mod_element 
        = XmlDomUtil::getChildElement(element, mod_element_name.c_str(), 0);
    mod_ptr_ = ModBase::getModPtrFromXml(mod_element);
  }
  std::string local_element_name = LocalAnno::getXmlElementName();;
  int local_count = XmlDomUtil::getChildCount(element, local_element_name.c_str());
  if (local_count != 0) {
    xercesc::DOMElement * local_element
        = XmlDomUtil::getChildElement(element, local_element_name.c_str(), 0); 
    local_anno_ptr_ = LocalAnnoPtr(new LocalAnno(local_element));
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
  if (mod_ptr_ != nullptr) {
    mod_ptr_->appendToXml(xml_doc, element);
  }
  if (local_anno_ptr_ != nullptr) {
    local_anno_ptr_->appendToXml(xml_doc, element);
  }
  parent->appendChild(element);
}

bool Change::cmpPosInc(const ChangePtr &a, const ChangePtr &b) {
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
  ModPtr mod_ptr = ori_ptr->mod_ptr_;
  ChangePtr change_ptr(
      new Change(left_bp_pos, right_bp_pos, change_type_ptr, mass_shift, mod_ptr));
  return change_ptr;
}

void Change::setLocalAnno(LocalAnnoPtr p) {
  local_anno_ptr_ = p;
  if (p != nullptr) {
    left_bp_pos_ = p->getLeftBpPos();
    right_bp_pos_ = p->getRightBpPos() + 1;
  }
}

}
