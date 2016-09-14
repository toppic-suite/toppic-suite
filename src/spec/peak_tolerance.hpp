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


#ifndef PROT_SPEC_PEAK_TOLERANCE_HPP_
#define PROT_SPEC_PEAK_TOLERANCE_HPP_

#include <memory>

#include "base/xml_dom_document.hpp"

namespace prot {

class PeakTolerance {
 public:
  PeakTolerance(double ppo, bool use_min_tolerance,
                double min_tolerance);

  PeakTolerance(xercesc::DOMElement* element);

  double compStrictErrorTole(double mass);

  // consider zero ptm relaxed error
  double compRelaxErrorTole(double m1, double m2) {
    return compStrictErrorTole(m1 + m2);
  }

  double getPpo() {return ppo_;}

  bool isUseMinTolerance() {return use_min_tolerance_;}

  double getMinTolerance() {return min_tolerance_;}

  void setPpo(double ppo) {ppo_ = ppo;}

  void setUseMinTolerance(bool use_min_tolerance) {
    use_min_tolerance_ = use_min_tolerance;
  }

  void setMinTolerance(double min_tolerance) {
    min_tolerance_ = min_tolerance;
  }

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "peak_tolerance";}

 private:
  double ppo_;
  /* whether or not use minimum tolerance */
  bool use_min_tolerance_;
  double min_tolerance_;
};

typedef std::shared_ptr<PeakTolerance> PeakTolerancePtr;

} /* namespace prot */

#endif    
