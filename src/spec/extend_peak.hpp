// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#ifndef PROT_EXTEND_PEAK_HPP_
#define PROT_EXTEND_PEAK_HPP_

#include <vector>
#include <memory>
#include <algorithm>

#include "spec/deconv_peak.hpp"

namespace prot {

class ExtendPeak;
typedef std::shared_ptr<ExtendPeak> ExtendPeakPtr;

class ExtendPeak : public Peak{
 public:
  ExtendPeak();

  ExtendPeak(DeconvPeakPtr base_peak_ptr, double mono_mass, double score);

  DeconvPeakPtr getBasePeakPtr(){return base_peak_ptr_;}

  double getMonoMass(){return mono_mass_;}

  double getScore(){return score_;}

  double getOrigTolerance(){return orig_tolerance_;}

  double getReverseTolerance(){return reverse_tolerance_;}

  void setOrigTolerance(double orig_tolerance) {
    orig_tolerance_ = orig_tolerance;}

  void setReverseTolerance(double reverse_tolerance) {
    reverse_tolerance_ = reverse_tolerance;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static bool cmpPosIncrease(const ExtendPeakPtr &a, const ExtendPeakPtr &b){
    return a->getPosition() < b->getPosition();
  }

  static std::string getXmlElementName() {return "extend_peak";}

 private:
  DeconvPeakPtr base_peak_ptr_;
  double mono_mass_;
  double score_;
  double orig_tolerance_;
  double reverse_tolerance_;
};

typedef std::vector<ExtendPeakPtr> ExtendPeakPtrVec;



} /* namespace prot */

#endif /* EXTEND_PEAK_HPP_ */
