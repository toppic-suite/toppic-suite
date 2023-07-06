//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <cmath>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/mass_constant.hpp"
#include "ms/spec/peak_util.hpp"
#include "ms/spec/precursor.hpp"

namespace toppic {

Precursor::Precursor(int id, double mono_mz, int charge, 
                     double inte, double apex_time):
    id_(id),
    mono_mz_(mono_mz),
    charge_(charge),
    inte_ (inte),
    apex_time_(apex_time) {}

Precursor::Precursor(XmlDOMElement* element) {
  id_ = xml_dom_util::getIntChildValue(element, "id", 0);
  mono_mz_ = xml_dom_util::getDoubleChildValue(element, "mono_mz", 0);
  charge_ = xml_dom_util::getIntChildValue(element, "charge", 0);
  inte_ = xml_dom_util::getDoubleChildValue(element, "intensity", 0);
  apex_time_ = xml_dom_util::getDoubleChildValue(element, "apex_time", 0);
}

double Precursor::getMonoMz() {
  if (std::isnan(mono_mz_)) {
    LOG_INFO("id " << id_ << " monoisotopic mz is not initialized!");
    return 0.0; 
  } else {
    return mono_mz_;
  }
}

double Precursor::getMonoMass() {
  if (mono_mz_ < 0 || std::isnan(mono_mz_)) {
    LOG_INFO("id " << id_ << " monoisotopic mass is not initialized!");
    return 0.0;
  } else {
    return peak_util::compPeakNeutralMass(mono_mz_, charge_);
  }
}

double Precursor::getMonoMassMinusWater() {
  if (mono_mz_ < 0 || std::isnan(mono_mz_)) {
    LOG_INFO("monoisotopic mass is not initialized!");
    return 0.0;
  } else {
    return peak_util::compPeakNeutralMass(mono_mz_, charge_)
        - mass_constant::getWaterMass();
  }
}

std::pair<int, int> Precursor::getMonoMassMinusWaterError(double ppo, double scale) {
  int mass = static_cast<int>(std::round(getMonoMassMinusWater() * scale));
  double error_tolerance = getErrorTolerance(ppo);
  int error = static_cast<int>(std::ceil(error_tolerance*scale));
  std::pair<int, int> result(mass, error);
  return result;
}

XmlDOMElement* Precursor::getPrecursorXml(XmlDOMDocument* xml_doc) {
  // float number precison
  int precison = 4;
  std::string precursor_str = Precursor::getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(precursor_str.c_str()); 
  std::string str = str_util::toString(id_);
  xml_doc->addElement(element, "id", str.c_str());
  str = str_util::fixedToString(mono_mz_, precison);
  xml_doc->addElement(element, "mono_mz", str.c_str());
  str = str_util::toString(charge_);
  xml_doc->addElement(element, "charge", str.c_str());
  str = str_util::toString(inte_);
  xml_doc->addElement(element, "intensity", str.c_str());
  str = str_util::toString(apex_time_);
  xml_doc->addElement(element, "apex_time", str.c_str());
  return element;
}

void Precursor::appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  parent->appendChild(getPrecursorXml(xml_doc));
}

bool Precursor::cmpInteDec(const PrecursorPtr &a, const PrecursorPtr &b) {
  return a->getInte() > b->getInte();
}

}  // namespace toppic
