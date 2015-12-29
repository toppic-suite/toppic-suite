#include "base/string_util.hpp"
#include "prsmview/anno_change.hpp"

namespace prot {

/*
AnnoExpectedChange::AnnoExpectedChange(ChangeTypePtr change_type_ptr, PtmPtr ptm_ptr) {
  change_type_ptr_ = change_type_ptr;
  ptm_ptr_ = ptm_ptr;
}
*/

void AnnoChange::addOccurence(int pos, const std::string &acid_letter) {
  std::pair<int, std::string> new_occurence(pos, acid_letter);
  occurences_.push_back(new_occurence);
}

void AnnoChange::appendExpectedXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("expected_change");
  std::string str = change_type_ptr_->getName();
  xml_doc->addElement(element, "change_type", str.c_str());
  PtmPtr ptm_ptr = mod_ptr_->getModResiduePtr()->getPtmPtr();
  ptm_ptr->appendAbbrNameToXml(xml_doc,element);

  for (size_t i = 0; i < occurences_.size(); i++) {
    xercesc::DOMElement* position_element = xml_doc->createElement("occurence");
    std::string str = StringUtil::convertToString(occurences_[i].first);
    xml_doc->addElement(position_element, "position", str.c_str());
    xml_doc->addElement(position_element, "acid_letter", occurences_[i].second.c_str());
    element->appendChild(position_element);
  }
  parent->appendChild(element);
}

void AnnoChange::appendUnexpectedXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent,
                                     int precison){
  xercesc::DOMElement* element = xml_doc->createElement("unexpected_change");
  std::string str = StringUtil::convertToString(left_bp_pos_);
  xml_doc->addElement(element, "left_position", str.c_str());
  str = StringUtil::convertToString(right_bp_pos_);
  xml_doc->addElement(element, "right_position", str.c_str());
  str = StringUtil::convertToString(mass_shift_, precison);
  xml_doc->addElement(element, "mass_shift", str.c_str());
  str = StringUtil::convertToString(color_);
  xml_doc->addElement(element, "unexpected_change_color", str.c_str());
  str = change_type_ptr_->getName();
  xml_doc->addElement(element, "type", str.c_str());
  
  /*
  if (type_ == "SHIFT") {
      xml_doc->addElement(element, "change_type", StringUtil::convertToString(4).c_str());
  }
  std::string occu = "";
  if (ptm_ptr_ != nullptr) {
      ptm_ptr_->appendAbbrNameToXml(xml_doc, element);
      for (size_t i = 0; i < occurences_.size(); i++) {
          if (i != occurences_.size() - 1) {
              occu += occurences_[i].second + StringUtil::convertToString(occurences_[i].first);
              occu += " / ";
          } else {
              occu += occurences_[i].second + StringUtil::convertToString(occurences_[i].first);
          }
      }
  } else {
      if (occurences_.size() > 0) {
          occu += occurences_[0].second + StringUtil::convertToString(occurences_[0].first);
          occu += " - ";
          occu += occurences_[1].second + StringUtil::convertToString(occurences_[1].first);
      }
  }
  xml_doc->addElement(element, "occurence", occu.c_str());
  */
  parent->appendChild(element);
}

AnnoChangePtr findExpectedChange(const AnnoChangePtrVec &expected_change_ptrs, 
                                 ChangeTypePtr change_type_ptr, ModPtr mod_ptr) {
  for (size_t i = 0; i < expected_change_ptrs.size(); i++) {
    if ((expected_change_ptrs[i]->getChangeTypePtr() == change_type_ptr) &&
        (expected_change_ptrs[i]->getModPtr() == mod_ptr)) {
      return expected_change_ptrs[i];
    }
  }
  return nullptr;
}

}
