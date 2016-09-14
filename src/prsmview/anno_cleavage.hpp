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


#ifndef PROT_ANNO_CLEAVAGE_HPP_
#define PROT_ANNO_CLEAVAGE_HPP_

#include "base/proteoform.hpp"
#include "spec/extend_peak.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "prsm/prsm.hpp"

namespace prot {

#define CLEAVAGE_TYPE_NORMAL "normal"
#define CLEAVAGE_TYPE_N_TRUNCATION "n_truncation"
#define CLEAVAGE_TYPE_C_TRUNCATION "c_truncation"
#define CLEAVAGE_TYPE_SEQ_START "seq_start"
#define CLEAVAGE_TYPE_SEQ_END "seq_end"

class AnnoCleavage {
 public:
  AnnoCleavage(int pos, const PeakIonPairPtrVec &pairs, bool exist_n_ion, bool exist_c_ion);

  void setPairs(PeakIonPairPtrVec pairs) {pairs_ = pairs;} 

  void setExistNIon(bool n) {exist_n_ion_ = n;};

  void setExistCIon(bool c) {exist_c_ion_ = c;};

  void setType(const std::string &type) {type_=type;}

  void setUnexpectedChange(bool u) {is_unexpected_change_ = u;}
  
  void setUnexpectedChangeColor(int color) {unexpected_change_color_ = color;}

  std::string getType(){return type_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  int pos_;
  bool exist_n_ion_;
  bool exist_c_ion_;
  std::string type_;
  PeakIonPairPtrVec pairs_;
  bool is_unexpected_change_;
  int unexpected_change_color_;
};

typedef std::shared_ptr<AnnoCleavage> AnnoCleavagePtr;
typedef std::vector<AnnoCleavagePtr> AnnoCleavagePtrVec;

AnnoCleavagePtrVec getProteoCleavage(PrsmPtr prsm_ptr, double min_mass);
} /* namespace prot */

#endif /* PROT_ANNO_CLEAVAGE_HPP_ */
