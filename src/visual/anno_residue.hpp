//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#ifndef TOPPIC_VISUAL_ANNO_RESIDUE_HPP_
#define TOPPIC_VISUAL_ANNO_RESIDUE_HPP_

#include "common/xml/xml_dom_document.hpp"
#include "common/base/residue.hpp"
#include "seq/proteoform.hpp"

namespace toppic {

class AnnoResidue;
typedef std::shared_ptr<AnnoResidue> AnnoResiduePtr;
typedef std::vector<AnnoResiduePtr> AnnoResiduePtrVec;

class AnnoResidue : public Residue {
 public:
  AnnoResidue(ResiduePtr residue_ptr, int pos);

  void appendViewXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static AnnoResiduePtrVec getAnnoResidues(ProteoformPtr proteoform_ptr);

 private:
  //residues position
  int pos_ = 0;
};

}  // namespace toppic

#endif /* TOPPIC_VISUAL_ANNO_RESIDUE_HPP_ */
