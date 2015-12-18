#include "base/string_util.hpp"
#include "prsmview/anno_expected_change.hpp"

namespace prot {

AnnoExpectedChange::AnnoExpectedChange(ChangeTypePtr change_type_ptr, PtmPtr ptm_ptr) {
  change_type_ptr_ = change_type_ptr;
  ptm_ptr_ = ptm_ptr;
}

void AnnoExpectedChange::addOccurence(int pos, const std::string &acid_letter) {
  std::pair<int, std::string> new_occurence(pos, acid_letter);
  occurences_.push_back(new_occurence);
}

void AnnoExpectedChange::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("expected_change");
  std::string str = change_type_ptr_->getName();
  xml_doc->addElement(element, "change_type", str.c_str());
  ptm_ptr_->appendAbbrNameToXml(xml_doc,element);

  for (size_t i = 0; i < occurences_.size(); i++) {
    xercesc::DOMElement* position_element = xml_doc->createElement("occurence");
    std::string str = StringUtil::convertToString(occurences_[i].first);
    xml_doc->addElement(position_element, "position", str.c_str());
    xml_doc->addElement(position_element, "acid_letter", occurences_[i].second.c_str());
    element->appendChild(position_element);
  }
  parent->appendChild(element);
}

AnnoExpectedChangePtr findExpectedChange(const AnnoExpectedChangePtrVec &expected_change_ptrs, 
                                         ChangeTypePtr change_type_ptr, PtmPtr ptm_ptr) {
  for (size_t i = 0; i < expected_change_ptrs.size(); i++) {
    if ((expected_change_ptrs[i]->getChangeTypePtr() == change_type_ptr) &&
        (expected_change_ptrs[i]->getPtmPtr() == ptm_ptr)) {
      return expected_change_ptrs[i];
    }
  }
  return nullptr;
}

}
