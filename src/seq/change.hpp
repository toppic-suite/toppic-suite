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


#ifndef TOPPIC_SEQ_CHANGE_HPP_
#define TOPPIC_SEQ_CHANGE_HPP_

#include "common/xml/xml_dom_element.hpp"
#include "common/base/mod.hpp"
#include "seq/mass_shift_type.hpp"
#include "seq/local_anno.hpp"

namespace toppic {

class XmlDOMDocument;

class Change;
typedef std::shared_ptr<Change> ChangePtr;

class Change {
 public:
  Change(int left_bp_pos, int right_bp_pos,
         MassShiftTypePtr type_ptr,
         double mass, ModPtr mod_ptr);

  explicit Change(XmlDOMElement* change_element);

  int getLeftBpPos() {return left_bp_pos_;}

  void setLeftBpPos(int p) {left_bp_pos_ = p;}

  int getRightBpPos() {return right_bp_pos_;}

  void setRightBpPos(int p) {right_bp_pos_ = p;}

  MassShiftTypePtr getTypePtr() {return type_ptr_;}

  double getMass() {return mass_;}

  void setMass(double m) {mass_ = m;}

  ModPtr getModPtr() {return mod_ptr_;}

  LocalAnnoPtr getLocalAnno() {return local_anno_ptr_;}

  void setLocalAnno(LocalAnnoPtr p);

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  static std::string getXmlElementName() {return "change";}

  static ChangePtr geneChangePtr(ChangePtr ori_change_ptr, int start_pos);

 protected:
  // left and right positions are based on break point positions
  int left_bp_pos_;

  int right_bp_pos_;

  MassShiftTypePtr type_ptr_;

  double mass_;

  ModPtr mod_ptr_;

  LocalAnnoPtr local_anno_ptr_;
};

typedef std::vector<ChangePtr> ChangePtrVec;

}  // namespace toppic

#endif

