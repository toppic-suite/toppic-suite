#include <sstream>
#include <stdlib.h>

#include <base/logger.hpp>
#include <base/string_util.hpp>
#include "spec/peak.hpp"
#include "spec/ms_header.hpp"

namespace prot {


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
  for (unsigned int i = 1; i < scans_.size(); i++) {
    scan_list <<  " " << scans_[i];
  }
  return scan_list.str();
}

void MsHeader::setScans(std::string s) {
  if (s == "") {
    scans_.clear();
    scans_.push_back(-1);
    return;
  }
  std::vector<std::string> strs;
  split(s,' ', strs); 
  for (unsigned int i = 0; i < strs.size(); i++) {
    scans_.push_back(atoi(strs[i].c_str()));
  }
}

void MsHeader::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
	xercesc::DOMElement* element = xml_doc->createElement("ms_header");
	xml_doc->addElement(element, "file_name", file_name_.c_str());
	std::string str = convertToString(id_);
	xml_doc->addElement(element, "id_", str.c_str());
	str = convertToString(prec_id_);
	xml_doc->addElement(element, "prec_id", str.c_str());
	xml_doc->addElement(element, "title", title_.c_str());
	str = convertToString(level_);
	xml_doc->addElement(element, "level", str.c_str());
	xercesc::DOMElement* scans = xml_doc->createElement("scan_list");
	for(unsigned int i=0;i< scans_.size();i++){
		str = convertToString(scans_[i]);
		xml_doc->addElement(scans, "scan", str.c_str());
	}
	str = convertToString(retention_time_);
	xml_doc->addElement(element, "retention_time", str.c_str());
	str = convertToString(prec_sp_mz_);
	xml_doc->addElement(element, "prec_sp_mz", str.c_str());
	str = convertToString(prec_mono_mz_);
	xml_doc->addElement(element, "prec_mono_mz_", str.c_str());
	str = convertToString(prec_charge_);
	xml_doc->addElement(element, "prec_charge_", str.c_str());
	str = convertToString(error_tolerance_);
	xml_doc->addElement(element, "error_tolerance_", str.c_str());
	activation_ptr_->appendXml(xml_doc,element);
	parent->appendChild(element);
}

}
