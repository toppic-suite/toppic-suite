#include <iostream>
#include <cmath>
#include <algorithm>

#include "base/logger.hpp"
#include "base/proteoform_factory.hpp"
#include "base/proteoform.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

#include "prsm/simple_prsm.hpp"

namespace prot {

SimplePrsm::SimplePrsm(MsHeaderPtr header_ptr, int spectrum_num, 
                       ProteoformPtr proteo_ptr,int score): 
  spectrum_num_(spectrum_num),
  //proteo_ptr_(proteo_ptr),
  score_(score) {
  spectrum_id_ = header_ptr->getId();
  spectrum_scan_ = header_ptr->getScansString();
  precursor_id_ = header_ptr->getPrecId();
  prec_mass_ = header_ptr->getPrecMonoMass();
  seq_name_ = proteo_ptr->getSeqName();
  seq_desc_ = proteo_ptr->getSeqDesc();
}

SimplePrsm::SimplePrsm(xercesc::DOMElement* element){
  spectrum_id_ = XmlDomUtil::getIntChildValue(element, "spectrum_id", 0);
  spectrum_scan_ = XmlDomUtil::getChildValue(element, "spectrum_scan", 0);
  precursor_id_ = XmlDomUtil::getIntChildValue(element, "precursor_id", 0);
  prec_mass_ = XmlDomUtil::getDoubleChildValue(element, "precursor_mass", 0);
  spectrum_num_ = XmlDomUtil::getDoubleChildValue(element, "spectrum_number", 0);
  seq_name_ = XmlDomUtil::getChildValue(element, "sequence_name", 0);
  seq_desc_ = XmlDomUtil::getChildValue(element, "sequence_desc", 0);
  score_ = XmlDomUtil::getDoubleChildValue(element, "score", 0);
}

xercesc::DOMElement* SimplePrsm::toXml(XmlDOMDocument* xml_doc){
  std::string element_name = SimplePrsm::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(spectrum_id_);
  xml_doc->addElement(element, "spectrum_id", str.c_str());
  xml_doc->addElement(element, "spectrum_scan", spectrum_scan_.c_str());
  str = StringUtil::convertToString(precursor_id_);
  xml_doc->addElement(element, "precursor_id", str.c_str());
  str = StringUtil::convertToString(prec_mass_);
  xml_doc->addElement(element, "precursor_mass", str.c_str());
  str = StringUtil::convertToString(spectrum_num_);
  xml_doc->addElement(element, "spectrum_number", str.c_str());
  str = StringUtil::convertToString(score_);
  xml_doc->addElement(element, "score", str.c_str());
  xml_doc->addElement(element, "sequence_name", seq_name_.c_str());
  xml_doc->addElement(element, "sequence_desc", seq_desc_.c_str());

  element_name = NTermShift::getXmlElementName() + "_list";
  xercesc::DOMElement* shift_list = xml_doc->createElement(element_name.c_str());
  for(size_t i=0; i<n_term_shifts_.size(); i++) {
    n_term_shifts_[i]->appendXml(xml_doc, shift_list);
  }
  element->appendChild(shift_list);
  return element;
}

} /* namespace prot */
