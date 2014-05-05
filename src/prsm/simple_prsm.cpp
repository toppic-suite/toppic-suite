/*
 * simple_prsm.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#include <iostream>
#include <cmath>

#include "simple_prsm.hpp"
#include "base/logger.hpp"
#include "base/proteoform.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

SimplePrsm::SimplePrsm(MsHeaderPtr header,ProteoformPtr seq,int score){
    spectrum_id_ = header->getId();
    spectrum_scan_ = header->getScansString();
    precursor_id_ = header->getPrecId();
    prec_mass_ = header->getPrecMonoMass();
    seq_= seq;
    seq_id_ = seq->getDbResSeqPtr()->getId();
    seq_name_ = seq->getDbResSeqPtr()->getName();
    score_ = score;
}

SimplePrsm::SimplePrsm(xercesc::DOMElement* element){
    spectrum_id_ = getIntChildValue(element, "spectrum_id", 0);
    spectrum_scan_ = getChildValue(element, "spectrum_scan", 0);
    precursor_id_ = getIntChildValue(element, "precursor_id", 0);
    prec_mass_ = getDoubleChildValue(element, "precursor_mass", 0);
    seq_id_ = getIntChildValue(element, "sequence_id", 0);
    seq_name_ = getChildValue(element, "sequence_name", 0);
    score_ = getDoubleChildValue(element, "score", 0);
}

int SimplePrsm::compareTo(SimplePrsmPtr simple_prsm_ptr){
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
void SimplePrsm::findSeq(std::vector<ProteoformPtr> seqs){
    seq_ = seqs[seq_id_];

    if(seq_->getSeqId() != seq_id_ || seq_->getSeqName() != seq_name_){
        std::cout<< "Sequence ID and/or name is not consistent!" << std::endl;
        std::exit(0);
    }
}

xercesc::DOMElement* SimplePrsm::toXml(XmlDOMDocument* xml_doc){
    xercesc::DOMElement* element = xml_doc->createElement("simple_prsm");
    std::string str = convertToString(spectrum_id_);
    xml_doc->addElement(element, "spectrum_id", str.c_str());
    xml_doc->addElement(element, "spectrum_scan", spectrum_scan_.c_str());
    str = convertToString(precursor_id_);
    xml_doc->addElement(element, "precursor_id", str.c_str());
    str = convertToString(prec_mass_);
    xml_doc->addElement(element, "precursor_mass", str.c_str());
    str = convertToString(seq_id_);
    xml_doc->addElement(element, "sequence_id", str.c_str());
    str = convertToString(score_);
    xml_doc->addElement(element, "score", str.c_str());
    xml_doc->addElement(element, "sequence_name", seq_name_.c_str());
    return element;
}
bool SimplePrsm::isMatch(MsHeaderPtr header){
    int header_spectrum_id = header->getId();
    std::string header_spectrum_scan = header->getScansString();
    int header_precursor_id = header->getPrecId();
    double header_precursor_mass = header->getPrecMonoMass();
    if(header_spectrum_id == spectrum_id_ 
     && header_precursor_id == precursor_id_){
        if(std::abs(header_precursor_mass-prec_mass_) > 0.00001 ||
                header_spectrum_scan != spectrum_scan_){
            LOG_ERROR("Error in combine simple Prsms! ");
        }
        return true;
    }
    else {
        return false;
    }
}

SimplePrsmPtrVec readSimplePrsm(std::string filename){
    SimplePrsmPtrVec results;
    XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
    if(parser){
        XmlDOMDocument* doc = new XmlDOMDocument(parser, filename.c_str());
        if (doc) {
            xercesc::DOMElement* root = doc->getDocumentElement();
            int simple_prsm_num = prot::getChildCount(root, "simple_prsm");
            for (int i = 0; i < simple_prsm_num; i++) {
                xercesc::DOMElement* simple_prsm 
            = getChildElement(root, "simple_prsm", i);
                results.push_back(SimplePrsmPtr(new SimplePrsm(simple_prsm)));
            }
        }
        delete doc;
    }
    return results;
}

SimplePrsmPtrVec findSimplePrsms(SimplePrsmPtrVec simple_prsm,
                                 MsHeaderPtr header){
    SimplePrsmPtrVec prsms ;
    for(unsigned int i=0;i<simple_prsm.size();i++){
        SimplePrsmPtr prsm = simple_prsm[i];
        if(prsm->isMatch(header)){
            prsms.push_back(prsm);
        }
    }
    return prsms;
}

} /* namespace prot */
