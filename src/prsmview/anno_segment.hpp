//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

