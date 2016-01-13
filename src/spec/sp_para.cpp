#include "base/activation_base.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "spec/sp_para.hpp"

namespace prot {

SpPara::SpPara(int min_peak_num,double min_mass,
               double min_extend_mass, 
               const std::vector<double> &ext_offsets,
               PeakTolerancePtr peak_tolerance_ptr,
               ActivationPtr activation_ptr): 
    min_peak_num_(min_peak_num),
    min_mass_(min_mass),
    extend_min_mass_(min_extend_mass),
    ext_offsets_(ext_offsets),
    peak_tolerance_ptr_(peak_tolerance_ptr),
    activation_ptr_(activation_ptr) {
    }

SpPara::SpPara(xercesc::DOMElement* element){
  min_peak_num_ = XmlDomUtil::getIntChildValue(element,"min_peak_num",0);
  min_mass_ = XmlDomUtil::getDoubleChildValue(element,"min_mass",0);
  extend_min_mass_ = XmlDomUtil::getDoubleChildValue(element,"extend_min_mass",0);
  xercesc::DOMElement* list_element 
      = XmlDomUtil::getChildElement(element, "extend_offset_list", 0);
  int offset_num =  XmlDomUtil::getChildCount(list_element, "extend_offset");
  for (int i = 0; i < offset_num; i++) {
    double offset = XmlDomUtil::getDoubleChildValue(list_element, "extend_offset", i);
    ext_offsets_.push_back(offset);
  }
  std::string element_name = PeakTolerance::getXmlElementName();
  xercesc::DOMElement* pt_element = XmlDomUtil::getChildElement(element,element_name.c_str(),0);
  peak_tolerance_ptr_ = PeakTolerancePtr(new PeakTolerance(pt_element));

  element_name = Activation::getXmlElementName();
  xercesc::DOMElement* ac_element = XmlDomUtil::getChildElement(element,element_name.c_str(),0);
  activation_ptr_ = ActivationBase::getActivationPtrFromXml(ac_element);
}

void SpPara::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = SpPara::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(min_peak_num_);
  xml_doc->addElement(element, "min_peak_num", str.c_str());
  str = StringUtil::convertToString(min_mass_);
  xml_doc->addElement(element, "min_mass", str.c_str());
  str = StringUtil::convertToString(extend_min_mass_).c_str();
  xml_doc->addElement(element, "extend_min_mass", str.c_str());
  xercesc::DOMElement* list_element 
      = xml_doc->createElement("extend_offset_list");
  element->appendChild(list_element);
  for (size_t i = 0; i < ext_offsets_.size(); i++) {
    str = StringUtil::convertToString(ext_offsets_[i]);
    xml_doc->addElement(list_element, "extend_offset", str.c_str());
  }
  element->appendChild(list_element);
  peak_tolerance_ptr_->appendXml(xml_doc, element);
  activation_ptr_->appendNameToXml(xml_doc,element);
  parent->appendChild(element); 
}

} /* namespace prot */
