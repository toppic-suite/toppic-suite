#include <base/logger.hpp>

#include "base/trunc.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

TruncPtrVec TruncFactory::trunc_ptr_vec_;

Trunc::Trunc(std::string name, int trunc_len, std::string str) {
  name_ = name;
  trunc_len_ = trunc_len;
  shift_ = 0;
  for (unsigned int i = 0; i < str.length(); i++) {
    std::string letter = str.substr(i, 1);
    AcidPtr acid_ptr = AcidFactory::getBaseAcidPtrByOneLetter(letter);
    acid_str_.push_back(acid_ptr);
    shift_ = shift_ - acid_ptr->getMonoMass();
  }
}

void Trunc::appendxml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
	xercesc::DOMElement* element = xml_doc->createElement("truncation");
	xml_doc->addElement(element, "name", name_.c_str());
	std::string str = convertToString(trunc_len_);
	xml_doc->addElement(element, "trunc_len", str.c_str());
	str = convertToString(shift_);
	xml_doc->addElement(element, "shift", str.c_str());
	xercesc::DOMElement* acid_list = xml_doc->createElement("amino_acid_list");
	for(unsigned int i=0;i<acid_str_.size();i++){
		acid_str_[i]->appendxml(xml_doc,acid_list);
	}
	element->appendChild(acid_list);
	parent->appendChild(element);
}

bool Trunc::isSameTrunc(int len,ResSeqPtr resseq){
	if(trunc_len_ != len){
		return false;
	}
	for(int i=0;i<trunc_len_;i++){
		if(acid_str_[i]->getName().compare(resseq->getResiduePtr(i)->getAcidPtr()->getName())!=0){
			return false;
		}
	}
	return true;
}



TruncPtr findProtTermTrunc(TruncPtrVec truncs,int trunc_len,ResSeqPtr resseq){
	for(unsigned int i=0;i<truncs.size();i++){
		if(truncs[i]->isSameTrunc(trunc_len,resseq)){
			return truncs[i];
		}
	}
	return nullptr;
};

TruncPtr findProtNTermTrunc(ResSeqPtr seq,int trunc_len,TruncPtrVec allowed_trunc){
//	ResSeqPtr resseq = seq->getResSeqPtr();
	return findProtTermTrunc(allowed_trunc,trunc_len,seq);
};
TruncPtr findProtCTermTrunc(ResSeqPtr seq,int last_res_pos,TruncPtrVec allowed_trunc){
//	int trunc_len = seq->getResSeqPtr()->getLen()-1-last_res_pos;
//	ResSeqPtr resseq = seq->getResSeqPtr();
	int trunc_len = seq->getLen()-1-last_res_pos;
	return  findProtTermTrunc(allowed_trunc,trunc_len,seq);
};

bool isAlignPrefix(TruncPtr n_trunc,double pep_n_term_shift,double threshold){
	if(n_trunc != nullptr && pep_n_term_shift <= threshold){
		return true;
	}
	return false;
};

bool isAlignSuffix(TruncPtr c_trunc,double pep_c_term_shift,double threshold){
	if(c_trunc != nullptr && pep_c_term_shift <= threshold){
		return true;
	}
	return false;
}

void TruncFactory::initFactory(const std::string &file_name) {
  prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int trunc_num = getChildCount(parent, "truncation");
    for (int i = 0; i < trunc_num; i++) {
      xercesc::DOMElement* element = getChildElement(parent, "truncation", i);
      std::string name = getChildValue(element, "name", 0);
      int trunc_len = getIntChildValue(element, "trunc_len", 0);
      std::string str = getChildValue(element, "acid_str", 0);
      LOG_DEBUG( "name " << name << " str " << str << " trunc len " << trunc_len);
      trunc_ptr_vec_.push_back(TruncPtr(new Trunc(name, trunc_len, str)));
    }
  }
}

TruncPtr TruncFactory::getBaseTruncPtrByName(const std::string &name) {
  for (unsigned int i = 0; i < trunc_ptr_vec_.size(); i++) {
    std::string n = trunc_ptr_vec_[i]->getName();
    if (n == name) {
      return trunc_ptr_vec_[i];
    }
  }
  return TruncPtr(nullptr);
}

}
