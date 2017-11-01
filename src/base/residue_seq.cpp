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


#include <sstream>

#include "base/residue_seq.hpp"
#include "base/string_util.hpp"

namespace prot {

ResidueSeq::ResidueSeq(const ResiduePtrVec &residues): 
  residues_(residues) {
  /* get residue mass sum */
  residue_mass_sum_ = 0;
  for (size_t i = 0; i < residues_.size(); i++) {
    residue_mass_sum_ += residues_[i]->getMass();
  }
}


ResSeqPtr ResidueSeq::getSubResidueSeq(int bgn, int end) {
  if (end - bgn < 0) {
    return getEmptyResidueSeq();
  } else {
    ResiduePtrVec sub_residues; 
    //from bgn to end,the sum of residues shoule be end - bgn + 1
    std::copy (residues_.begin() + bgn, residues_.begin() + end + 1,
               std::back_inserter(sub_residues) );
    return ResSeqPtr(new ResidueSeq(sub_residues));
  }
}

std::string ResidueSeq::toString() {
  std::stringstream s;
  for (size_t i = 0; i < residues_.size(); i++) {
    s << residues_[i]->toString();
  }
  s<< std::endl;
  return s.str();
}

std::string ResidueSeq::toAcidString() {
  std::stringstream s;
  for (size_t i = 0; i < residues_.size(); i++) {
    s << residues_[i]->getAcidPtr()->getOneLetter();
  }
  return s.str();
}

/*
void ResidueSeq::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = ResidueSeq::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(residue_mass_sum_);
  xml_doc->addElement(element, "residue_mass_sum", str.c_str());
  xercesc::DOMElement* residuelist = xml_doc->createElement("residue_list");
  for(size_t i=0;i<residues_.size();i++){
    residues_[i]->appendXml(xml_doc,residuelist);
  }
  element->appendChild(residuelist);
  parent->appendChild(element);
}
*/

ResSeqPtr ResidueSeq::getEmptyResidueSeq() {
  ResiduePtrVec residues;
  return ResSeqPtr(new ResidueSeq(residues));
}

}
