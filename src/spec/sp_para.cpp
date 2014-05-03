/*
 * sp_para.cpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#include <spec/sp_para.hpp>

namespace prot {

SpPara::SpPara(int min_peak_num,double min_mass,
               const PeakTolerancePtr &peak_tolerance,
               const ExtendSpParaPtr &extend_sp_para,
               const ActivationPtr &activation){
  activation_ = activation;
  min_peak_num_ = min_peak_num;
  min_mass_ = min_mass;
  peak_tolerance_ptr_=peak_tolerance;
  extend_sp_para_ = extend_sp_para;
}

SpPara::SpPara(xercesc::DOMElement* element){
  min_peak_num_ = prot::getIntChildValue(element,"min_peak_num",0);
  min_mass_ = prot::getDoubleChildValue(element,"min_mass",0);
  extend_sp_para_ = ExtendSpParaPtr(new ExtendSpPara(getChildElement(element,"extend_sp_para",0)));
  peak_tolerance_ptr_ = PeakTolerancePtr(new PeakTolerance(getChildElement(element,"peak_tolerance",0)));
}

void SpPara::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = xml_doc->createElement("sp_para");
  xml_doc->addElement(element, "min_peak_num", prot::convertToString(min_peak_num_).c_str());
  xml_doc->addElement(element, "min_mass", prot::convertToString(min_mass_).c_str());
  extend_sp_para_->appendXml(xml_doc, element);
  peak_tolerance_ptr_->appendXml(xml_doc, element);
  parent->appendChild(element); 
}
} /* namespace prot */
