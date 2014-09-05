
#include "base/anno_residue.hpp"

namespace prot {

AnnoResidue::AnnoResidue(ResiduePtr residue_ptr, int display_pos):
    Residue(residue_ptr->getAcidPtr(),residue_ptr->getPtmPtr()) {
      display_pos_= display_pos;
      type_= ANNO_RESIDUE_TYPE_NORMAL;
      is_unexpected_change_ = false;
      unexpected_change_color_ = 0;
    }

void AnnoResidue::appendViewXml(XmlDOMDocument* xml_doc,
                                xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("character");
  xml_doc->addElement(element, "type", "residue");
  std::string str = convertToString(display_pos_);
  xml_doc->addElement(element, "display_position", str.c_str());
  str = getAcidPtr()->getOneLetter();
  xml_doc->addElement(element, "acid", str.c_str());
  str = type_;
  xml_doc->addElement(element, "residue_type", str.c_str());
  str = convertToString(is_unexpected_change_);
  xml_doc->addElement(element, "is_unexpected_change", str.c_str());
  str = convertToString(unexpected_change_color_);
  xml_doc->addElement(element, "unexpected_change_color_", str.c_str());
  parent->appendChild(element);
}

}

