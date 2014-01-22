#include <sstream>

#include "base/residue_seq.hpp"

namespace prot {

ResidueSeq::ResidueSeq(ResiduePtrVec residues) {
  residues_ = residues;
  /* get residue mass sum */
  residue_mass_sum_ = 0;
  for (unsigned int i = 0; i < residues_.size(); i++) {
    residue_mass_sum_ += residues_[i]->getMass();
  }
}


ResSeqPtr ResidueSeq::getSubResidueSeq(int bgn, int end) {
  if (end - bgn < 0) {
    return getEmptyResidueSeq();
  } else {
    ResiduePtrVec sub_residues; 
    std::copy (residues_.begin() + bgn, residues_.begin() + end, 
               std::back_inserter(sub_residues) );
    return ResSeqPtr( new ResidueSeq(sub_residues));
  }
}

std::string ResidueSeq::toString() {
  std::stringstream s;
  for (unsigned int i = 0; i < residues_.size(); i++) {
    s << residues_[i]->toString();
  }
  s<< std::endl;
  return s.str();
}

void ResidueSeq::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
	xercesc::DOMElement* element = xml_doc->createElement("residue_seq");
	std::string str = convertToString(residue_mass_sum_);
	xml_doc->addElement(element, "residue_mass_sum", str.c_str());
	xercesc::DOMElement* residuelist = xml_doc->createElement("residue_list");
	for(unsigned int i=0;i<residues_.size();i++){
		residues_[i]->appendXml(xml_doc,residuelist);
	}
	element->appendChild(residuelist);
	parent->appendChild(element);
}

ResSeqPtr getEmptyResidueSeq() {
  ResiduePtrVec residues;
  return ResSeqPtr(new ResidueSeq(residues));
}

}
