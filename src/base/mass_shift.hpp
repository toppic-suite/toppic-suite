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

#ifndef PROT_BASE_MASS_SHIFT_HPP_
#define PROT_BASE_MASS_SHIFT_HPP_

#include <memory>
#include <vector>
#include <string>

#include "base/change.hpp"

namespace prot {

class MassShift;
typedef std::shared_ptr<MassShift> MassShiftPtr;

class MassShift {
 public:
  MassShift(int left_bp_pos, int right_bp_pos, MassShiftTypePtr type_ptr):
      left_bp_pos_(left_bp_pos),
      right_bp_pos_(right_bp_pos),
      type_ptr_(type_ptr),
      shift_(0.0) { }

  explicit MassShift(xercesc::DOMElement* mass_shift_element);

  int getLeftBpPos() {return left_bp_pos_;}

  void setLeftBpPos(int p) {left_bp_pos_ = p;}

  int getRightBpPos() {return right_bp_pos_;}

  void setRightBpPos(int p) {right_bp_pos_ = p;}

  std::vector<ChangePtr> getChangePtrVec() {return change_vec_;}

  ChangePtr getChangePtr(size_t idx) {return change_vec_[idx];}

  void setChangePtr(ChangePtr change);

  MassShiftTypePtr getTypePtr() {return type_ptr_;}

  double getMassShift() {return shift_;}

  std::string getSeqStr();

  static std::string getXmlElementName() {return "mass_shift";}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static bool cmpPosInc(const MassShiftPtr & a, const MassShiftPtr & b);

 private:
  int left_bp_pos_;

  int right_bp_pos_;

  MassShiftTypePtr type_ptr_;

  double shift_;

  std::vector<ChangePtr> change_vec_;
};

typedef std::vector<MassShiftPtr> MassShiftPtrVec;

}  // namespace prot

#endif
