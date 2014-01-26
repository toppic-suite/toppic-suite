
#include "spec/extend_sp_para.hpp"

namespace prot {

ExtendSpPara::ExtendSpPara(double extend_min_mass, 
                           std::vector<double> ext_offsets) {
  extend_min_mass_ = extend_min_mass;
  ext_offsets_ = ext_offsets;
}

ExtendSpPara::ExtendSpPara(xercesc::DOMElement* element) {
  extend_min_mass_ = getDoubleChildValue(element, "extend_min_mass", 0);
  xercesc::DOMElement* list_element 
      = getChildElement(element, "extend_offset_list", 0);
  int offset_num =  getChildCount(list_element, "extend_offset");
  for (int i = 0; i < offset_num; i++) {
    double offset = getDoubleChildValue(list_element, "extend_offset", i);
    ext_offsets_.push_back(offset);
  }
}

void ExtendSpPara::appendXml(XmlDOMDocument* xml_doc, 
                             xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = xml_doc->createElement("extend_sp_para");
  std::string str = convertToString(extend_min_mass_);
  xml_doc->addElement(element, "extend_min_mass", str.c_str());
  xercesc::DOMElement* list_element 
      = xml_doc->createElement("extend_offset_list");
  element->appendChild(list_element);
  for (unsigned int i = 0; i < ext_offsets_.size(); i++) {
    str = convertToString(ext_offsets_[i]);
    xml_doc->addElement(list_element, "extend_offset", str.c_str());
  }
  parent->appendChild(element);
}

}
