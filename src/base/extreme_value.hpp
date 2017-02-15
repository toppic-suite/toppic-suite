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


#ifndef PROT_BASE_EXTREME_VALUE_HPP_
#define PROT_BASE_EXTREME_VALUE_HPP_

#include <memory>
#include <vector>
#include "base/xml_dom_document.hpp"

namespace prot {

class ExtremeValue;
typedef std::shared_ptr<ExtremeValue> ExtremeValuePtr;

class ExtremeValue {
 public:
  ExtremeValue (double one_prot_prob, double test_num, double adjust_factor);

  ExtremeValue (xercesc::DOMElement* element);

  double getPValue() {return p_value_;}

  double getEValue() {return e_value_;}

  double getOneProtProb() {return one_prot_prob_;}

  double getTestNum() {return test_num_;}

  double getAdjustFactor() { return adjust_factor_;}

  void setOneProtProb(double one_prot_prob);

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "extreme_value";}

  static ExtremeValuePtr getMaxEvaluePtr();

 private:
  // one_prot_prob is the probability that the spectrum and a randem problem 
  // have a protein-spectrum-match with a score no less than the threshold 
  double one_prot_prob_;
  double test_num_;
  double adjust_factor_;
  double p_value_;
  double e_value_;

  void init();
};

typedef std::vector<ExtremeValuePtr> ExtremeValuePtrVec;

}

#endif
