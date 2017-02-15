// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "base/logger.hpp"
#include "base/string_util.hpp"
#include "base/algorithm.hpp"
#include "spec/theo_peak.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

PeakIonPair::PeakIonPair(MsHeaderPtr ms_header_ptr, ExtendPeakPtr real_peak_ptr, 
                         TheoPeakPtr theo_peak_ptr): 
    ms_header_ptr_(ms_header_ptr),
    real_peak_ptr_(real_peak_ptr),
    theo_peak_ptr_(theo_peak_ptr) {
    }

void PeakIonPair::appendRealPeakToXml(XmlDOMDocument* xml_doc, 
                                      xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = xml_doc->createElement("matched_peak");
  std::string str = theo_peak_ptr_->getIonPtr()->getIonTypePtr()->getName();
  xml_doc->addElement(element, "ion_type", str.c_str());
  str = StringUtil::convertToString(theo_peak_ptr_->getIonPtr()->getPos());
  xml_doc->addElement(element, "ion_position", str.c_str());
  str = StringUtil::convertToString(theo_peak_ptr_->getIonPtr()->getDisplayPos());
  xml_doc->addElement(element, "ion_display_position", str.c_str());
  str = StringUtil::convertToString(ms_header_ptr_->getId());
  xml_doc->addElement(element, "spec_id", str.c_str());
  str = StringUtil::convertToString(real_peak_ptr_->getBasePeakPtr()->getId());
  xml_doc->addElement(element, "peak_id", str.c_str());
  str = StringUtil::convertToString(real_peak_ptr_->getBasePeakPtr()->getCharge());
  xml_doc->addElement(element, "peak_charge", str.c_str());
  parent->appendChild(element);
}

void PeakIonPair::appendTheoPeakToXml(XmlDOMDocument* xml_doc, 
                                      xercesc::DOMElement* parent) {
  int precison=4;
  xercesc::DOMElement* element = xml_doc->createElement("matched_ion");
  std::string str 
      = theo_peak_ptr_->getIonPtr()->getIonTypePtr()->getName();
  xml_doc->addElement(element, "ion_type", str.c_str());
  str = StringUtil::convertToString(theo_peak_ptr_->getShift());
  xml_doc->addElement(element, "match_shift", str.c_str()); 
  str = StringUtil::convertToString(theo_peak_ptr_->getModMass(), precison);
  xml_doc->addElement(element, "theoretical_mass", str.c_str()); 
  str = StringUtil::convertToString(theo_peak_ptr_->getIonPtr()->getPos());
  xml_doc->addElement(element, "ion_position", str.c_str());
  str = StringUtil::convertToString(theo_peak_ptr_->getIonPtr()->getDisplayPos());
  xml_doc->addElement(element, "ion_display_position", str.c_str());
  str = theo_peak_ptr_->getIonPtr()->getIonTypePtr()->getName();
  // convert display position to a string with five letters.
  std::string disp_pos = StringUtil::convertToString(theo_peak_ptr_->getIonPtr()->getDisplayPos());
  while (disp_pos.length() < 5) {
    disp_pos = "0" + disp_pos;
  }
  str += disp_pos;
  xml_doc->addElement(element, "ion_sort_name", str.c_str());
  str = StringUtil::convertToString(theo_peak_ptr_->getIonPtr()->getPos());
  xml_doc->addElement(element, "ion_left_position", str.c_str());
  double error = real_peak_ptr_->getMonoMass() - theo_peak_ptr_->getModMass();
  str = StringUtil::convertToString(error,precison);
  xml_doc->addElement(element, "mass_error", str.c_str()); 
  str = StringUtil::convertToString(error * 1000000 / real_peak_ptr_->getMonoMass(),precison-2);
  xml_doc->addElement(element, "ppm", str.c_str()); 
  parent->appendChild(element);
}

bool PeakIonPair::cmpRealPeakPosInc(const PeakIonPairPtr &a, 
                                    const PeakIonPairPtr &b) {
  return a->getRealPeakPtr()->getBasePeakPtr()->getPosition() 
      < b->getRealPeakPtr()->getBasePeakPtr()->getPosition();

}

}
