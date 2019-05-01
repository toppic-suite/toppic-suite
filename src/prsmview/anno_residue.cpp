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


#include <string>

#include "common/base/residue_util.hpp"
#include "common/util/str_util.hpp"
#include "prsmview/anno_residue.hpp"

namespace toppic {

AnnoResidue::AnnoResidue(ResiduePtr residue_ptr, int pos):
    Residue(residue_ptr->getAminoAcidPtr(), residue_ptr->getPtmPtr()),
    pos_(pos) {}

void AnnoResidue::appendViewXml(XmlDOMDocument* xml_doc,
                                xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = xml_doc->createElement("residue");
  std::string str = str_util::toString(pos_);
  xml_doc->addElement(element, "position", str.c_str());

  str = getAminoAcidPtr()->getOneLetter();
  xml_doc->addElement(element, "acid", str.c_str());

  parent->appendChild(element);
}

AnnoResiduePtrVec AnnoResidue::getAnnoResidues(ProteoformPtr proteoform_ptr) {
  StringPairVec acid_ptm_pairs 
      = proteoform_ptr->getFastaSeqPtr()->getAcidPtmPairVec();
  ResiduePtrVec fasta_residues 
      = residue_util::convertStrToResiduePtrVec(acid_ptm_pairs);

  AnnoResiduePtrVec res_ptrs;
  int prot_len = fasta_residues.size();
  for (int i = 0; i < prot_len; i++) {
    res_ptrs.push_back(std::make_shared<AnnoResidue>(fasta_residues[i], i));
  }
  return res_ptrs;
}


}  // namespace toppic

