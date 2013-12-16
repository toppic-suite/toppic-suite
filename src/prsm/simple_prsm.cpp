/*
 * simple_prsm.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#include <log4cxx/logger.h>
#include <iostream>

#include "simple_prsm.hpp"
#include "base/proteoform.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

log4cxx::LoggerPtr sPrsm_logger(log4cxx::Logger::getLogger("SimplePrSM"));

SimplePrSM::SimplePrSM(MsHeaderPtr header,ProteoformPtr seq,int score){
	spectrum_id_ = header->getId();
	spectrum_scan_ = header->getScansString();
	precursor_id_ = header->getPrecId();
	prec_mass_ = header->getPrecMonoMass();
	seq_= seq;
	//todo:seq_id and name should be delete? name
//	seq_id_ = seq->getResSeqPtr()->getResiduePtr();
//	seq_name_ = seq->getResSeqPtr();
	score_ = score;
}

SimplePrSM::SimplePrSM(xercesc::DOMElement* element){
	spectrum_id_ = getIntChildValue(element, "spectrum_id", 0);
	spectrum_scan_ = getChildValue(element, "spectrum_scan", 0);
	precursor_id_ = getIntChildValue(element, "precursor_id", 0);
	prec_mass_ = getDoubleChildValue(element, "precursor_mass", 0);
//	seq_= getIntChildValue(element, "spectrum_id_", 0);
	seq_id_ = getIntChildValue(element, "sequence_id", 0);
	seq_name_ = getChildValue(element, "sequence_name", 0);
	score_ = getDoubleChildValue(element, "score", 0);
}

int SimplePrSM::compareTo(SimplePrSMPtr simple_prsm_ptr){
	if (simple_prsm_ptr->getScore()>getScore()){
		return 1;
	}
	else if (simple_prsm_ptr->getScore() < getScore()){
		return -1;
	}
	else {
		return 0;
	}
}
void SimplePrSM::findSeq(std::vector<ProteoformPtr> seqs){
//todo: didn't find id
}

xercesc::DOMElement* SimplePrSM::toXml(XmlDOMDocument* xml_doc){
	xercesc::DOMElement* element = xml_doc->createElement("simple_prsm");
	std::string str = convertToString(spectrum_id_);
	xml_doc->addElement(element, "spectrum_id", str.c_str());
	xml_doc->addElement(element, "spectrum_scan", spectrum_scan_.c_str());
	str = convertToString(precursor_id_);
	xml_doc->addElement(element, "precursor_id", str.c_str());
	str = convertToString(prec_mass_);
	xml_doc->addElement(element, "prec_mass", str.c_str());
	str = convertToString(seq_id_);
	xml_doc->addElement(element, "seq_id", str.c_str());
	str = convertToString(score_);
	xml_doc->addElement(element, "score", str.c_str());
	xml_doc->addElement(element, "seq_name", seq_name_.c_str());
	return element;
}
bool SimplePrSM::isMatch(MsHeaderPtr header){
	int header_spectrum_id = header->getId();
	std::string header_spectrum_scan = header->getScansString();
	int header_precursor_id = header->getPrecId();
	double header_precursor_mass = header->getPrecMonoMass();
	if(header_spectrum_id == spectrum_id_ && header_precursor_id == precursor_id_){
		if(header_precursor_mass != prec_mass_ ||
				header_spectrum_scan.compare(spectrum_scan_)!=0){
			LOG4CXX_ERROR(sPrsm_logger, "Error in combine simple PrSMs! ");
		}
		return true;
	}
	else {
		return false;
	}
}

SimplePrSMPtrVec readSimplePrSM(const char * filename){
	SimplePrSMPtrVec results;
	XmlDOMParser* parser = getXmlDOMInstance();
	if(parser){
		std::cout<< filename << std::endl;
		XmlDOMDocument* doc = new XmlDOMDocument(parser, filename);
		if (doc) {
			xercesc::DOMElement* root = doc->getDocumentElement();
			xercesc::DOMElement* parent = getChildElement(root, "simple_prsm_list", 0);
			int simple_prsm_num = prot::getChildCount(parent, "simple_prsm");
			for (int i = 0; i < simple_prsm_num; i++) {
				xercesc::DOMElement* simple_prsm = getChildElement(parent, "simple_prsm", i);
				results.push_back(SimplePrSMPtr(new SimplePrSM(simple_prsm)));
			}
		}
		delete doc;
	}
	return results;
}

} /* namespace prot */
