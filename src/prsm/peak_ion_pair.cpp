//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include <string>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace toppic {

PeakIonPair::PeakIonPair(MsHeaderPtr ms_header_ptr, ExtendPeakPtr real_peak_ptr,
                         TheoPeakPtr theo_peak_ptr): 
    ms_header_ptr_(ms_header_ptr),
    real_peak_ptr_(real_peak_ptr),
    theo_peak_ptr_(theo_peak_ptr) {}

void PeakIonPair::appendRealPeakToXml(XmlDOMDocument* xml_doc, 
                                      XmlDOMElement* parent) {
  XmlDOMElement* element = xml_doc->createElement("matched_peak");
  std::string str = theo_peak_ptr_->getIonPtr()->getIonTypePtr()->getName();
  xml_doc->addElement(element, "ion_type", str.c_str());
  str = str_util::toString(theo_peak_ptr_->getIonPtr()->getPos());
  xml_doc->addElement(element, "ion_position", str.c_str());
  str = str_util::toString(theo_peak_ptr_->getIonPtr()->getDisplayPos());
  xml_doc->addElement(element, "ion_display_position", str.c_str());
  str = str_util::toString(ms_header_ptr_->getId());
  xml_doc->addElement(element, "spec_id", str.c_str());
  str = str_util::toString(real_peak_ptr_->getBasePeakPtr()->getId());
  xml_doc->addElement(element, "peak_id", str.c_str());
  str = str_util::toString(real_peak_ptr_->getBasePeakPtr()->getCharge());
  xml_doc->addElement(element, "peak_charge", str.c_str());
  parent->appendChild(element);
}

void PeakIonPair::appendTheoPeakToXml(XmlDOMDocument* xml_doc, 
                                      XmlDOMElement* parent) {
  int precison = 4;
  XmlDOMElement* element = xml_doc->createElement("matched_ion");
  std::string str 
      = theo_peak_ptr_->getIonPtr()->getIonTypePtr()->getName();
  xml_doc->addElement(element, "ion_type", str.c_str());
  str = str_util::toString(theo_peak_ptr_->getShift());
  xml_doc->addElement(element, "match_shift", str.c_str()); 
  str = str_util::fixedToString(theo_peak_ptr_->getModMass(), precison);
  xml_doc->addElement(element, "theoretical_mass", str.c_str()); 
  str = str_util::toString(theo_peak_ptr_->getIonPtr()->getPos());
  xml_doc->addElement(element, "ion_position", str.c_str());
  str = str_util::toString(theo_peak_ptr_->getIonPtr()->getDisplayPos());
  xml_doc->addElement(element, "ion_display_position", str.c_str());
  str = theo_peak_ptr_->getIonPtr()->getIonTypePtr()->getName();
  // convert display position to a string with five letters.
  std::string disp_pos = str_util::toString(theo_peak_ptr_->getIonPtr()->getDisplayPos());
  while (disp_pos.length() < 5) {
    disp_pos = "0" + disp_pos;
  }
  str += disp_pos;
  xml_doc->addElement(element, "ion_sort_name", str.c_str());
  str = str_util::toString(theo_peak_ptr_->getIonPtr()->getPos());
  xml_doc->addElement(element, "ion_left_position", str.c_str());
  double error = real_peak_ptr_->getMonoMass() - theo_peak_ptr_->getModMass();
  str = str_util::fixedToString(error, precison);
  xml_doc->addElement(element, "mass_error", str.c_str()); 
  str = str_util::fixedToString(error * 1000000 / real_peak_ptr_->getMonoMass(), precison - 2);
  xml_doc->addElement(element, "ppm", str.c_str()); 
  parent->appendChild(element);
}

bool PeakIonPair::cmpRealPeakPosInc(const PeakIonPairPtr &a, const PeakIonPairPtr &b) {
  return a->getRealPeakPtr()->getBasePeakPtr()->getPosition() 
      < b->getRealPeakPtr()->getBasePeakPtr()->getPosition();
}

bool PeakIonPair::cmpTheoPeakPosInc(const PeakIonPairPtr &a, const PeakIonPairPtr &b) {
  return a->getTheoPeakPtr()->getIonPtr()->getPos() 
      < b->getTheoPeakPtr()->getIonPtr()->getPos();
}

}  // namespace toppic
