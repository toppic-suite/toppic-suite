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


#ifndef PROT_ANNO_SEGMENT_HPP_
#define PROT_ANNO_SEGMENT_HPP_

#include "base/change.hpp"
#include "base/change_type.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class AnnoSegment {
 public:
  AnnoSegment(std::string segment_type, int left_pos, int right_pos,
              double mass_shift, int color);

  std::string getType() {return segment_type_;}

  int getRightPos() {return right_pos_;}

  void setRightPos(int right_pos) {right_pos_ = right_pos;}

  void addOccurence(int pos, const std::string &acid_letter);

  void setPtmPtr(PtmPtr p) {ptm_ptr_ = p;}

  void addScr(double s) {score_.push_back(s);}

  std::string getResidueAnno();

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent, int precison);

 private:
  std::string segment_type_;
  std::string anno_;
  std::string occu_;
  int left_pos_;
  int right_pos_;
  double mass_shift_;
  int color_;
  PtmPtr ptm_ptr_;
  std::vector<std::pair<int, std::string>> occurences_;
  std::vector<double> score_;
};

typedef std::shared_ptr<AnnoSegment> AnnoSegmentPtr;
typedef std::vector<AnnoSegmentPtr> AnnoSegmentPtrVec;

}
#endif

