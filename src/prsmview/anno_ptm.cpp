#include "base/string_util.hpp"
#include "prsmview/anno_ptm.hpp"

namespace prot {

AnnoPtm::AnnoPtm(PtmPtr ptm_ptr, ChangeTypePtr change_type_ptr) {
  ptm_ptr_ = ptm_ptr;
  change_type_ptr_ = change_type_ptr;
}

void AnnoPtm::addOccurence(int pos, const std::string &acid_letter) {
  std::pair<int, std::string> new_occurence(pos, acid_letter);
  occurences_.push_back(new_occurence);
}

void AnnoPtm::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
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

AnnoPtmPtr AnnoPtm::findPtm(const AnnoPtmPtrVec &ptm_ptrs, PtmPtr ptm_ptr, 
                            ChangeTypePtr change_type_ptr) {
  for (size_t i = 0; i < ptm_ptrs.size(); i++) {
    if ((ptm_ptrs[i]->getChangeTypePtr() == change_type_ptr) &&
        (ptm_ptrs[i]->getPtmPtr() == ptm_ptr)) {
      return ptm_ptrs[i];
    }
  }
  return nullptr;
}

}
