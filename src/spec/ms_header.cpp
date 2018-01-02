//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#include <sstream>
#include <cmath>
#include <utility>
#include <string>
#include <vector>

#include "base/activation_base.hpp"
#include "base/logger.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "spec/peak_util.hpp"
#include "spec/ms_header.hpp"

namespace prot {

MsHeader::MsHeader(xercesc::DOMElement* element) {
  file_name_ = xml_dom_util::getChildValue(element, "file_name", 0);
  id_ = xml_dom_util::getIntChildValue(element, "id", 0);
  prec_id_ = xml_dom_util::getIntChildValue(element, "prec_id", 0);
  title_ = xml_dom_util::getChildValue(element, "title", 0);
  level_ = xml_dom_util::getIntChildValue(element, "level", 0);

  xercesc::DOMElement* scan_element
      = xml_dom_util::getChildElement(element, "scan_list", 0);
  int scans = xml_dom_util::getChildCount(scan_element, "scan");
  for (int i = 0; i < scans; i++) {
    scans_.push_back(xml_dom_util::getIntChildValue(scan_element, "scan", i));
  }
  retention_time_ = xml_dom_util::getDoubleChildValue(element, "retention_time", 0);
  prec_sp_mz_ = xml_dom_util::getDoubleChildValue(element, "prec_sp_mz", 0);
  prec_mono_mz_ = xml_dom_util::getDoubleChildValue(element, "prec_mono_mz", 0);
  prec_charge_ = xml_dom_util::getIntChildValue(element, "prec_charge", 0);
  std::string element_name = Activation::getXmlElementName();
  xercesc::DOMElement* ac_element
      = xml_dom_util::getChildElement(element, element_name.c_str(), 0);
  activation_ptr_ = ActivationBase::getActivationPtrFromXml(ac_element);
}

double MsHeader::getPrecMonoMass() {
  if (prec_mono_mz_ < 0 || std::isnan(prec_mono_mz_)) {
    LOG_WARN("monoisotopic mass is not initialized");
    return 0.0;
  } else {
    return peak_util::compPeakMass(prec_mono_mz_, prec_charge_);
  }
}

double MsHeader::getPrecSpMass() {
  if (prec_sp_mz_ < 0) {
    LOG_WARN("precursor spectrum mass is not initialized");
    return 0.0;
  } else {
    return peak_util::compPeakMass(prec_sp_mz_, prec_charge_);
  }
}

double MsHeader::getPrecMonoMassMinusWater() {
  if (prec_mono_mz_ < 0 || std::isnan(prec_mono_mz_)) {
    LOG_WARN("monoisotopic mass is not initialized");
    return 0.0;
  } else {
    return peak_util::compPeakMass(prec_mono_mz_, prec_charge_)
        - mass_constant::getWaterMass();
  }
}

std::pair<int, int> MsHeader::getPrecMonoMassMinusWaterError(double ppo, double scale) {
  int mass = static_cast<int>(std::round(getPrecMonoMassMinusWater() * scale));
  double error_tolerance = getErrorTolerance(ppo);
  int error = static_cast<int>(std::ceil(error_tolerance*scale));
  std::pair<int, int> result(mass, error);
  return result;
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
  std::vector<std::string> strs = string_util::split(s, ' ');
  for (size_t i = 0; i < strs.size(); i++) {
    scans_.push_back(std::stoi(strs[i]));
  }
}

xercesc::DOMElement* MsHeader::getHeaderXml(XmlDOMDocument* xml_doc) {
  // float number precison
  int precison = 4;
  xercesc::DOMElement* element = xml_doc->createElement("ms_header");
  xml_doc->addElement(element, "file_name", file_name_.c_str());
  std::string str = string_util::convertToString(id_);
  xml_doc->addElement(element, "id", str.c_str());
  str = string_util::convertToString(prec_id_);
  xml_doc->addElement(element, "prec_id", str.c_str());
  xml_doc->addElement(element, "title", title_.c_str());
  str = string_util::convertToString(level_);
  xml_doc->addElement(element, "level", str.c_str());
  str = getScansString();
  xml_doc->addElement(element, "scans", str.c_str());
  xercesc::DOMElement* scans = xml_doc->createElement("scan_list");
  for (size_t i = 0; i < scans_.size(); i++) {
    str = string_util::convertToString(scans_[i]);
    xml_doc->addElement(scans, "scan", str.c_str());
  }
  element->appendChild(scans);
  str = string_util::convertToString(retention_time_);
  xml_doc->addElement(element, "retention_time", str.c_str());
  str = string_util::convertToString(prec_sp_mz_);
  xml_doc->addElement(element, "prec_sp_mz", str.c_str());
  str = string_util::convertToString(prec_mono_mz_, precison);
  xml_doc->addElement(element, "prec_mono_mz", str.c_str());
  str = string_util::convertToString(prec_charge_);
  xml_doc->addElement(element, "prec_charge", str.c_str());
  activation_ptr_->appendNameToXml(xml_doc, element);
  return element;
}

void MsHeader::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  parent->appendChild(getHeaderXml(xml_doc));
}

MsHeaderPtr MsHeader::geneMsHeaderPtr(MsHeaderPtr ori_ptr, double new_prec_mass) {
  MsHeaderPtr new_header_ptr = std::make_shared<MsHeader>(*ori_ptr.get());
  double mono_mz = peak_util::compMonoMz(new_prec_mass, ori_ptr->getPrecCharge());
  new_header_ptr->setPrecMonoMz(mono_mz);
  return new_header_ptr;
}

bool MsHeader::cmpPrecInteDec(const MsHeaderPtr &a, const MsHeaderPtr &b) {
  return a->getPrecInte() > b->getPrecInte();
}

}  // namespace prot
