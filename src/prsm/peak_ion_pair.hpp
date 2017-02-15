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


#ifndef PROT_PRSM_PEAK_ION_PAIR_HPP_
#define PROT_PRSM_PEAK_ION_PAIR_HPP_

#include "base/xml_dom_document.hpp"
#include "spec/ms_header.hpp"
#include "spec/extend_peak.hpp"
#include "spec/theo_peak.hpp"

namespace prot {

class PeakIonPair;
typedef std::shared_ptr<PeakIonPair> PeakIonPairPtr;

class PeakIonPair {
 public:
  PeakIonPair(MsHeaderPtr ms_header_ptr, ExtendPeakPtr real_peak_ptr, 
              TheoPeakPtr theo_peak_ptr); 

  MsHeaderPtr getMsHeaderPtr() {return ms_header_ptr_;}

  ExtendPeakPtr getRealPeakPtr() {return real_peak_ptr_;}

  TheoPeakPtr getTheoPeakPtr() {return theo_peak_ptr_;}

  void setId(int id) {id_ = id;}

  void appendRealPeakToXml(XmlDOMDocument* xml_doc, 
                           xercesc::DOMElement* parent);

  void appendTheoPeakToXml(XmlDOMDocument* xml_doc, 
                           xercesc::DOMElement* parent);

  static bool cmpRealPeakPosInc(const PeakIonPairPtr &a, const PeakIonPairPtr &b);

 private:
  int id_;
  MsHeaderPtr ms_header_ptr_;
  ExtendPeakPtr real_peak_ptr_;
  TheoPeakPtr theo_peak_ptr_;
};
typedef std::vector<PeakIonPairPtr> PeakIonPairPtrVec;
typedef std::vector<PeakIonPairPtrVec> PeakIonPairPtrVec2D;

}

#endif
