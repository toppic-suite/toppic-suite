/*
 * simple_prsm.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#include <log4cxx/logger.h>

#include "simple_prsm.hpp"
#include "proteoform.hpp"
#include "xml_dom_document.hpp"

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
//todo:should the element be shared_prt? peak_tolarance
xercesc::DOMElement* toXml(){
	//todo:implement toXml
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

} /* namespace prot */
