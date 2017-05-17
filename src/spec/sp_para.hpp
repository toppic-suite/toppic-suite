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


#ifndef PROT_SPEC_SP_PARA_HPP_
#define PROT_SPEC_SP_PARA_HPP_

#include <memory>
#include "base/activation.hpp"
#include "base/xml_dom_document.hpp"
#include "spec/peak_tolerance.hpp"

namespace prot {

class SpPara {
 public:
  SpPara(int min_peak_num, double min_mass_, double min_extend_mass, 
         const std::vector<double> &ext_offsets,
         PeakTolerancePtr peak_tolerance_ptr,
         ActivationPtr activation_ptr);

  SpPara(xercesc::DOMElement* element);

  double getMinMass() {return min_mass_;}

  double getExtendMinMass() {return extend_min_mass_;}

  const std::vector<double>& getExtendOffsets() {return ext_offsets_;}

  PeakTolerancePtr getPeakTolerancePtr() {return peak_tolerance_ptr_;}

  void setPeakTolerancePtr(PeakTolerancePtr peak_tolerance_ptr){
    peak_tolerance_ptr_ = peak_tolerance_ptr;}

  ActivationPtr getActivationPtr() {return activation_ptr_;}

  void setActivationPtr(ActivationPtr activation_ptr) {
    activation_ptr_ = activation_ptr;}

  int getMinPeakNum() {return min_peak_num_;}

  void setMinPeakNum(int min_peak_num) {min_peak_num_=min_peak_num;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "sp_para";}

  std::vector<double> mod_mass_;

 private:
  int min_peak_num_;

  // if the mass if smaller than min_mass, the mass is removed. 
  double min_mass_;
  // if the mass is smaller than extend_min_mass, the peak is not extended 
  double extend_min_mass_;
  std::vector<double> ext_offsets_;

  PeakTolerancePtr peak_tolerance_ptr_;
  ActivationPtr activation_ptr_;

};

typedef std::shared_ptr<SpPara> SpParaPtr;

} /* namespace prot */

#endif /* SP_PARA_HPP_ */
