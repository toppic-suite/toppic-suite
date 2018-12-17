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


#ifndef PROT_PRSM_PEAK_ION_PAIR_HPP_
#define PROT_PRSM_PEAK_ION_PAIR_HPP_

#include <vector>

#include "base/xml_dom_document.hpp"
#include "spec/ms_header.hpp"
#include "spec/extend_peak.hpp"
#include "spec/theo_peak.hpp"

namespace toppic {

class PeakIonPair;
typedef std::shared_ptr<PeakIonPair> PeakIonPairPtr;

class PeakIonPair {
 public:
  PeakIonPair(MsHeaderPtr ms_header_ptr, ExtendPeakPtr real_peak_ptr,
              TheoPeakPtr theo_peak_ptr): 
      ms_header_ptr_(ms_header_ptr),
      real_peak_ptr_(real_peak_ptr),
      theo_peak_ptr_(theo_peak_ptr) {}

  MsHeaderPtr getMsHeaderPtr() {return ms_header_ptr_;}

  ExtendPeakPtr getRealPeakPtr() {return real_peak_ptr_;}

  TheoPeakPtr getTheoPeakPtr() {return theo_peak_ptr_;}

  void setId(int id) {id_ = id;}

  void appendRealPeakToXml(XmlDOMDocument* xml_doc, 
                           xercesc::DOMElement* parent);

  void appendTheoPeakToXml(XmlDOMDocument* xml_doc, 
                           xercesc::DOMElement* parent);

  static bool cmpRealPeakPosInc(const PeakIonPairPtr &a, const PeakIonPairPtr &b);

  static bool cmpTheoPeakPosInc(const PeakIonPairPtr &a, const PeakIonPairPtr &b);

 private:
  int id_;

  MsHeaderPtr ms_header_ptr_;

  ExtendPeakPtr real_peak_ptr_;

  TheoPeakPtr theo_peak_ptr_;
};
typedef std::vector<PeakIonPairPtr> PeakIonPairPtrVec;
typedef std::vector<PeakIonPairPtrVec> PeakIonPairPtrVec2D;

}  // namespace toppic

#endif
