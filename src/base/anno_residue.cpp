/*
 * anno_residue.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: xunlikun
 */

#include "base/anno_residue.hpp"
namespace prot {

AnnoResidue::AnnoResidue(ResiduePtr residue_ptr):
    Residue(residue_ptr->getAcidPtr(),residue_ptr->getPtmPtr()) {
      pos_=0;
      type_="normal";
      display_pos_=0;
      is_modified_ = false;
      shift_=0;
      is_expected_ = false;
      display_bg_ = -1;
    }

void AnnoResidue::appendViewXml(XmlDOMDocument* xml_doc,
                                xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("character");
  xml_doc->addElement(element, "type", "residue");
  std::string str = convertToString(pos_);
  xml_doc->addElement(element, "position", str.c_str());
  str = getAcidPtr()->getOneLetter();
  xml_doc->addElement(element, "acid", str.c_str());
  str = type_;
  xml_doc->addElement(element, "residue_type", str.c_str());
  str = shift_style_;
  xml_doc->addElement(element, "shift_style", str.c_str());
  str = convertToString(is_expected_);
  xml_doc->addElement(element, "is_expected", str.c_str());
  str = convertToString(is_modified_);
  xml_doc->addElement(element, "is_modification", str.c_str());
  str = convertToString(shift_,2);
  xml_doc->addElement(element, "shift", str.c_str());
  str = convertToString(display_pos_);
  xml_doc->addElement(element, "display_position", str.c_str());
  str = convertToString(display_bg_);
  xml_doc->addElement(element, "display_background", str.c_str());
  parent->appendChild(element);
}

}

