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


ResidueSeq ResidueSeq::getSubResidueSeq(int bgn, int end) {
  if (end - bgn < 0) {
    return getEmptyResidueSeq();
  } else {
    ResiduePtrVec sub_residues; 
    std::copy (residues_.begin() + bgn, residues_.begin() + end, 
               std::back_inserter(sub_residues) );
    return ResidueSeq(sub_residues);
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
	for(int i=0;i<residues_.size();i++){
		residues_[i]->appendXml(xml_doc,residuelist);
	}
	element->appendChild(residuelist);
	parent->appendChild(element);
}

bool ResidueSeq::allowsMod(ProtModPtr mod){
	if(mod->getName().compare("NONE")==0){
		return true;
	}
	else if(mod->getName().compare("NME")==0){
		if(residues_.size()>=2 && residues_[0]->getAcidPtr()->getOneLetter().compare("M")==0){
			return true;
		}
		return false;
	}
	else if(mod->getName().compare("ACETYLATION")==0){
		return true;
	}
	else if(mod->getName().compare("NME_ACETYLATION")==0){
		if(residues_.size()>=2 && residues_[0]->getAcidPtr()->getOneLetter().compare("M")==0){
			return true;
		}
		return false;
	}

}

ResidueSeq getEmptyResidueSeq() {
  ResiduePtrVec residues;
  return ResidueSeq(residues);
}

}
