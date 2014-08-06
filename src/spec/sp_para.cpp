#include "spec/sp_para.hpp"

namespace prot {

SpPara::SpPara(int min_peak_num,double min_mass,
               double min_extend_mass, 
               const std::vector<double> &ext_offsets,
               PeakTolerancePtr peak_tolerance_ptr,
               ActivationPtr activation_ptr){
  min_peak_num_ = min_peak_num;
  min_mass_ = min_mass;
  extend_min_mass_ = min_extend_mass;
  ext_offsets_ = ext_offsets;
  peak_tolerance_ptr_=peak_tolerance_ptr;
  activation_ptr_ = activation_ptr;
}

SpPara::SpPara(xercesc::DOMElement* element){
  min_peak_num_ = getIntChildValue(element,"min_peak_num",0);
  min_mass_ = getDoubleChildValue(element,"min_mass",0);
  extend_min_mass_ = getDoubleChildValue(element,"extend_min_mass",0);
  xercesc::DOMElement* list_element 
      = getChildElement(element, "extend_offset_list", 0);
  int offset_num =  getChildCount(list_element, "extend_offset");
  for (int i = 0; i < offset_num; i++) {
    double offset = getDoubleChildValue(list_element, "extend_offset", i);
    ext_offsets_.push_back(offset);
  }
  peak_tolerance_ptr_ = PeakTolerancePtr(new PeakTolerance(getChildElement(element,"peak_tolerance",0)));
  std::string activation_name = getChildValue(element, "activation",0);
  activation_ptr_ = ActivationFactory::getBaseActivationPtrByName(activation_name);
}

void SpPara::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = xml_doc->createElement("sp_para");
  xml_doc->addElement(element, "min_peak_num", prot::convertToString(min_peak_num_).c_str());
  xml_doc->addElement(element, "min_mass", prot::convertToString(min_mass_).c_str());
  xml_doc->addElement(element, "extend_min_mass", prot::convertToString(extend_min_mass_).c_str());
  xercesc::DOMElement* list_element 
      = xml_doc->createElement("extend_offset_list");
  element->appendChild(list_element);
  for (size_t i = 0; i < ext_offsets_.size(); i++) {
    std::string str = convertToString(ext_offsets_[i]);
    xml_doc->addElement(list_element, "extend_offset", str.c_str());
  }
  element->appendChild(list_element);
  peak_tolerance_ptr_->appendXml(xml_doc, element);
  xml_doc->addElement(element, "activation", activation_ptr_->getName().c_str());
  parent->appendChild(element); 
}

} /* namespace prot */
