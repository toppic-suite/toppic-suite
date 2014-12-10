#include <sstream>
#include <stdlib.h>

#include <base/logger.hpp>
#include <base/string_util.hpp>
#include "spec/peak.hpp"
#include "spec/ms_header.hpp"

namespace prot {

MsHeader::MsHeader(xercesc::DOMElement* element){
    file_name_ = getChildValue(element,"file_name",0);
    id_ = getIntChildValue(element,"id",0);
    prec_id_ = getIntChildValue(element,"prec_id",0);
    title_= getChildValue(element,"title",0);
    level_ = getIntChildValue(element,"level",0);

    xercesc::DOMElement* scan_element 
      = prot::getChildElement(element,"scan_list",0);
    int scans = getChildCount(scan_element,"scan");
    for(int i=0;i<scans;i++){
        scans_.push_back(getIntChildValue(scan_element,"scan",i));
    }
    retention_time_ = getDoubleChildValue(element,"retention_time",0);
    prec_sp_mz_ = getDoubleChildValue(element,"prec_sp_mz",0);
    prec_mono_mz_ = getDoubleChildValue(element,"prec_mono_mz",0);
    prec_charge_ = getIntChildValue(element,"prec_charge",0);
    error_tolerance_ = getDoubleChildValue(element,"error_tolerance",0);
    xercesc::DOMElement* act_element 
      = prot::getChildElement(element,"activation",0);
    std::string act_name = getChildValue(act_element,"name",0);
    activation_ptr_= ActivationFactory::getBaseActivationPtrByName(act_name);
}

double MsHeader::getPrecMonoMass() {
  if (prec_mono_mz_ < 0) {
    LOG_WARN("monoisotopic mass is not initialized")
    return 0.0;
  } 
  else {
    return compPeakMass(prec_mono_mz_, prec_charge_);
  }
}

double MsHeader::getPrecSpMass() {
  if (prec_sp_mz_ < 0) {
    LOG_WARN("precursor spectrum mass is not initialized");
    return 0.0;
  } else {
    return compPeakMass(prec_sp_mz_, prec_charge_);
  }
}

double MsHeader::getPrecMonoMassMinusWater() {
  if (prec_mono_mz_ < 0) {
    LOG_WARN("monoisotopic mass is not initialized");
    return 0.0;
  } else {
    return compPeakMass(prec_mono_mz_, prec_charge_)
        - MassConstant::getWaterMass();
  }
}

std::string MsHeader::toString() {
  std::stringstream tmp;
        tmp << "MS HEADER\n";
        tmp << "==========\n";
        tmp << "Title = " << title_ << "\n";
        tmp << "Scan Number = " << scans_[0] << "\n";
        tmp << "MS Level = " << level_ << "\n";
        tmp << "Activation type = " << activation_ptr_ << "\n";
        tmp << "Precursor Sp Mz = " << prec_sp_mz_ << "\n";
        tmp << "Precursor Charge = " << prec_charge_ << "\n";
        tmp << "Precursro Mono Mz = " << prec_mono_mz_ << "\n";
        return tmp.str();
    }

std::string MsHeader::getScansString() {
  std::stringstream scan_list;
  scan_list << scans_[0];
  for (size_t i = 1; i < scans_.size(); i++) {
    scan_list <<  " " << scans_[i];
  }
  return scan_list.str();
}

void MsHeader::setScans(const std::string &s) {
  if (s == "") {
    scans_.clear();
    scans_.push_back(-1);
    return;
  }
  std::vector<std::string> strs;
  split(s,' ', strs); 
  for (size_t i = 0; i < strs.size(); i++) {
    scans_.push_back(std::stoi(strs[i]));
  }
}

xercesc::DOMElement* MsHeader::getHeaderXml(XmlDOMDocument* xml_doc) {
  int pos = 4;
  xercesc::DOMElement* element = xml_doc->createElement("ms_header");
  xml_doc->addElement(element, "file_name", file_name_.c_str());
  std::string str = convertToString(id_);
  xml_doc->addElement(element, "id", str.c_str());
  str = convertToString(prec_id_);
  xml_doc->addElement(element, "prec_id", str.c_str());
  xml_doc->addElement(element, "title", title_.c_str());
  str = convertToString(level_);
  xml_doc->addElement(element, "level", str.c_str());
  str = getScansString();
  xml_doc->addElement(element, "scans", str.c_str());
  xercesc::DOMElement* scans = xml_doc->createElement("scan_list");
  for (size_t i = 0; i < scans_.size(); i++) {
    str = convertToString(scans_[i]);
    xml_doc->addElement(scans, "scan", str.c_str());
  }
  element->appendChild(scans);
  str = convertToString(retention_time_);
  xml_doc->addElement(element, "retention_time", str.c_str());
  str = convertToString(prec_sp_mz_);
  xml_doc->addElement(element, "prec_sp_mz", str.c_str());
  str = convertToString(prec_mono_mz_);
  xml_doc->addElement(element, "prec_mono_mz", str.c_str());
  str = convertToString(prec_charge_);
  xml_doc->addElement(element, "prec_charge", str.c_str());
  str = convertToString(error_tolerance_);
  xml_doc->addElement(element, "error_tolerance", str.c_str());
  activation_ptr_->appendXml(xml_doc, element);
  str = convertToString(prec_mono_mz_, pos);
  xml_doc->addElement(element, "precursor_mz", str.c_str());
  if (prec_mono_mz_ < 0) {
    std::cout << "monoisotopic mass is not initialized!" << std::endl;
  } else {
    double prec_mass = prec_mono_mz_ * prec_charge_
        - prec_charge_ * prot::MassConstant::getProtonMass();
    str = convertToString(prec_mass, pos);
    xml_doc->addElement(element, "precursor_mass", str.c_str());
  }
  return element;
}

void MsHeader::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
    parent->appendChild(getHeaderXml(xml_doc));
}

}
