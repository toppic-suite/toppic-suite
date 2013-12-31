#include "base/change.hpp"

namespace prot {

Change::Change(int left_bp_pos, int right_bp_pos, int change_type,
         double mass_shift, PtmPtr ptm_ptr) {
  left_bp_pos_ = left_bp_pos;
  right_bp_pos_ = right_bp_pos;
  change_type_ = change_type;
  mass_shift_ = mass_shift;
  ptm_ptr_ = ptm_ptr;
}

Change::Change(Change ori, int start) {
  left_bp_pos_ = ori.left_bp_pos_ - start;
  right_bp_pos_ = ori.right_bp_pos_ - start;
  change_type_ = ori.change_type_;
  mass_shift_ = ori.mass_shift_;
  ptm_ptr_ = ori.ptm_ptr_;
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
	ptm_ptr_->appendxml(xml_doc,element);
	parent->appendChild(element);
}

}

