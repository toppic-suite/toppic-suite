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

#ifndef TOPPIC_VISUAL_ANNO_PTM_POSITION_HPP_
#define TOPPIC_VISUAL_ANNO_PTM_POSITION_HPP_

#include <memory>
#include <string>
#include <vector>

#include "common/xml/xml_dom_document.hpp"

namespace toppic {

class AnnoPtmPosition {
 public:
  AnnoPtmPosition(int left_pos, int right_pos, std::string anno);

  int getLeftPos() {return left_pos_;}

  int getRightPos() {return right_pos_;}

  std::string getAnno() {return anno_;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

 private:
  int left_pos_;
  
  int right_pos_;

  std::string anno_;
};

typedef std::shared_ptr<AnnoPtmPosition> AnnoPtmPositionPtr;
typedef std::vector<AnnoPtmPositionPtr> AnnoPtmPositionPtrVec;

}  // namespace toppic

#endif

