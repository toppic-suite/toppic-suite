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


#ifndef PROT_ANNO_RESIDUE_HPP_
#define PROT_ANNO_RESIDUE_HPP_

#include "base/residue.hpp"


namespace prot {

#define ANNO_RESIDUE_TYPE_NORMAL "normal"
#define ANNO_RESIDUE_TYPE_N_TRUNCATION "n_truncation"
#define ANNO_RESIDUE_TYPE_C_TRUNCATION "c_truncation"
#define ANNO_RESIDUE_TYPE_KNOWN_CHANGE "known_change"
#define ANNO_RESIDUE_TYPE_UNKNOWN_CHANGE "unknown_change"

class AnnoResidue : public Residue {
 public:
  AnnoResidue(ResiduePtr residue_ptr, int pos);

  void setType(const std::string &type) {type_ = type;}

  void setUnexpectedChange(bool u) {is_unexpected_change_ = u;}

  void setUnexpectedChangeColor(int color) {unexpected_change_color_ = color;}

  void setPossiblePosColor(int c) {possible_pos_color_ = c;}

  void setAnno(std::string s) {anno_ = s;}

  void appendViewXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

 private:
  //residues position
  int pos_ = 0;
  //residues type
  std::string type_ = ANNO_RESIDUE_TYPE_NORMAL;
  //is expected
  bool is_unexpected_change_ = false;
  //unexpected change color
  int unexpected_change_color_ = 0;
  //possible postion for known ptms
  int possible_pos_color_ = 0;
  // all possible positions with score
  std::string anno_;
};

typedef std::shared_ptr<AnnoResidue> AnnoResiduePtr;
typedef std::vector<AnnoResiduePtr> AnnoResiduePtrVec;

}

#endif /* PROT_ANNO_RESIDUE_HPP_ */
