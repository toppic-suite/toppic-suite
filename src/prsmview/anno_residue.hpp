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


#ifndef TOPPIC_PRSM_VIEW_ANNO_RESIDUE_HPP_
#define TOPPIC_PRSM_VIEW_ANNO_RESIDUE_HPP_

#include <string>
#include <vector>

#include "common/xml/xml_dom_document.hpp"
#include "common/base/residue.hpp"

namespace toppic {

const std::string ANNO_RESIDUE_TYPE_NORMAL         = "normal";

const std::string ANNO_RESIDUE_TYPE_N_TRUNCATION   = "n_truncation";

const std::string ANNO_RESIDUE_TYPE_C_TRUNCATION   = "c_truncation";

const std::string ANNO_RESIDUE_TYPE_KNOWN_CHANGE   = "known_change";

const std::string ANNO_RESIDUE_TYPE_UNKNOWN_CHANGE = "unknown_change";

class AnnoResidue : public Residue {
 public:
  AnnoResidue(ResiduePtr residue_ptr, int pos):
      Residue(residue_ptr->getAminoAcidPtr(), residue_ptr->getPtmPtr()),
      pos_(pos),
      type_(ANNO_RESIDUE_TYPE_NORMAL),
      is_unexpected_change_(false),
      unexpected_change_color_(0) {}

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

}  // namespace toppic

#endif /* TOPPIC_PRSM_VIEW_ANNO_RESIDUE_HPP_ */
