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

#ifndef TOPPIC_VISUAL_ANNO_MASS_SHIFT_HPP_
#define TOPPIC_VISUAL_ANNO_MASS_SHIFT_HPP_

#include <string>
#include <vector>

#include "common/xml/xml_dom_document.hpp"
#include "seq/alter_type.hpp"

namespace toppic {

class AnnoMassShift {
 public:
  AnnoMassShift(int id, int left_pos, int right_pos, 
                const std::string & anno_str, 
                AlterTypePtr & mass_shift_type);

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

 private:
  int id_;

  int left_pos_;

  int right_pos_;

  std::string anno_str_;

  AlterTypePtr mass_shift_type_;
};

typedef std::shared_ptr<AnnoMassShift> AnnoMassShiftPtr;

typedef std::vector<AnnoMassShiftPtr> AnnoMassShiftPtrVec;

}  // namespace toppic

#endif

