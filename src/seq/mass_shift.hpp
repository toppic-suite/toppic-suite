//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_SEQ_MASS_SHIFT_HPP_
#define TOPPIC_SEQ_MASS_SHIFT_HPP_

#include "seq/alter.hpp"

namespace toppic {

class MassShift;
typedef std::shared_ptr<MassShift> MassShiftPtr;

class MassShift {
 public:
  explicit MassShift(int left_bp_pos, int right_bp_pos, double shift);

  explicit MassShift(AlterPtr alter_ptr);  

  explicit MassShift(MassShiftPtr shift_ptr, int start);

  explicit MassShift(XmlDOMElement* mass_shift_element);

  int getLeftBpPos() {return left_bp_pos_;}

  void setLeftBpPos(int p) {left_bp_pos_ = p;}

  int getRightBpPos() {return right_bp_pos_;}

  void setRightBpPos(int p) {right_bp_pos_ = p;}

  std::vector<AlterPtr> getAlterPtrVec() {return alter_vec_;}

  AlterPtr getAlterPtr(size_t idx) {return alter_vec_[idx];}

  int getAlterNum() {return alter_vec_.size();}

  void setAlterPtrVec(AlterPtrVec alter_vec) {alter_vec_ = alter_vec;}

  AlterTypePtr getTypePtr(); 

  double getMassShift() {return shift_;}

  void setMassShift(double shift) {shift_ = shift;}

  std::string getAnnoStr();

  static std::string getXmlElementName() {return "mass_shift";}

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  static bool cmpPosInc(const MassShiftPtr & a, const MassShiftPtr & b);

 private:
  int left_bp_pos_;

  int right_bp_pos_;

  double shift_;

  // A mass shift is annotated by a list of alterations. 
  std::vector<AlterPtr> alter_vec_;
};

typedef std::vector<MassShiftPtr> MassShiftPtrVec;

}  // namespace toppic

#endif
