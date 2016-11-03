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


#ifndef PROT_BASE_LOCAL_ANNO_HPP_
#define PROT_BASE_LOCAL_ANNO_HPP_

#include <vector>
#include <algorithm>
#include "base/ptm.hpp"
#include "base/ptm_base.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class LocalAnno {
 public:
  LocalAnno(xercesc::DOMElement* element);

  LocalAnno(int left_pos, int right_pos, double conf, 
            std::vector<double> scr_vec, double raw_scr, PtmPtr p):
      left_pos_(left_pos), right_pos_(right_pos), 
      conf_(conf), scr_vec_(scr_vec), raw_scr_(raw_scr), ptm_ptr_(p){}

  int getLeftBpPos() {return left_pos_;}

  int getRightBpPos() {return right_pos_;}
  // get the confidence score
  double getConf() {return conf_;}
  // get the score vector containing score for each site
  std::vector<double> getScrVec() {return scr_vec_;}

  double getScr() {return std::accumulate(scr_vec_.begin(), scr_vec_.end(), 0.0);}

  PtmPtr getPtmPtr() {return ptm_ptr_;}

  double getRawScr() {return raw_scr_;};

  void getRawScr(double s) {raw_scr_ = s;};

  double getMassShift() {return mass_shift_;}

  void setMassShift(double m) {mass_shift_ = m;}

  static std::string getXmlElementName() {return "localization_annotation";}

  void appendToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  // left and right position filtered by thredshould
  int left_pos_, right_pos_;
  // confidence score
  double conf_;
  // score for each possible site
  std::vector<double> scr_vec_;
  // raw score
  double raw_scr_;
  // PTM
  PtmPtr ptm_ptr_;
  // mass shift
  double mass_shift_;
};

typedef std::shared_ptr<LocalAnno> LocalAnnoPtr;
typedef std::vector<LocalAnnoPtr> LocalAnnoPtrVec;

} // namespace prot

#endif
