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


#ifndef PROT_SPEC_DECONV_PEAK_HPP_
#define PROT_SPEC_DECONV_PEAK_HPP_

#include <memory>
#include <vector>
#include "spec/peak.hpp"

namespace prot {

class DeconvPeak;
typedef std::shared_ptr<DeconvPeak> DeconvPeakPtr;

class DeconvPeak : public Peak {
 public:
  DeconvPeak (int id, double mono_mass, double intensity, int charge);

  DeconvPeak(xercesc::DOMElement* element);

  int getCharge() {return charge_;}

  int getId() {return id_;}

  double getMonoMass() {return getPosition();}

  double getMonoMz() {return compMonoMz(getMonoMass(), charge_);}

  double getScore() {return score_;}

  void setId(int id) {id_ = id;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static bool cmpPosIncreasep(const DeconvPeakPtr &a, const DeconvPeakPtr &b){
    return a->getPosition() < b->getPosition();
  }

  static std::string getXmlElementName() {return "deconv_peak";}

 private:
  int id_;
  int charge_;
  double score_ = 1.0;
};

typedef std::vector<DeconvPeakPtr> DeconvPeakPtrVec;

}
#endif
