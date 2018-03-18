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


#ifndef PROT_ANNO_SEGMENT_HPP_
#define PROT_ANNO_SEGMENT_HPP_


#include <utility>
#include <string>
#include <vector>

#include "base/mass_shift.hpp"
#include "base/mass_shift_type.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class AnnoSegment {
 public:
  AnnoSegment(std::string segment_type, int left_pos, int right_pos,
              const std::string & match_seq, int color):
      segment_type_(segment_type),
      left_pos_(left_pos),
      right_pos_(right_pos),
      match_seq_(match_seq),
      color_(color),
      mass_shift_type_(MassShiftType::UNEXPECTED) {}

  std::string getType() {return segment_type_;}

  void setMassShiftType(MassShiftTypePtr shift_type) {
    mass_shift_type_ = shift_type;
  }

  int getRightPos() {return right_pos_;}

  void setRightPos(int right_pos) {right_pos_ = right_pos;}

  void addOccurence(int pos, const std::string &acid_letter);

  void setPtmPtr(PtmPtr p) {ptm_ptr_ = p;}

  void addScr(double s) {score_.push_back(s);}

  std::string getResidueAnno();

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

 private:
  // EMPTY or SHIFT
  std::string segment_type_;

  std::string anno_;

  int left_pos_;

  int right_pos_;

  std::string match_seq_;

  int color_;

  PtmPtr ptm_ptr_;

  std::vector<std::pair<int, std::string> > occurences_;

  std::vector<double> score_;

  MassShiftTypePtr mass_shift_type_;
};

typedef std::shared_ptr<AnnoSegment> AnnoSegmentPtr;

typedef std::vector<AnnoSegmentPtr> AnnoSegmentPtrVec;

}  // namespace prot

#endif

