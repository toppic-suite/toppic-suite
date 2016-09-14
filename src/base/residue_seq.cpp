// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
