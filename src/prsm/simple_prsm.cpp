#include <iostream>
#include <cmath>
#include <algorithm>

#include "simple_prsm.hpp"
#include "base/logger.hpp"
#include "base/proteoform_reader.hpp"
#include "base/proteoform.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

SimplePrsm::SimplePrsm(MsHeaderPtr header_ptr, int spectrum_num, 
                       ProteoformPtr proteo_ptr,int score){
  spectrum_id_ = header_ptr->getId();
  spectrum_scan_ = header_ptr->getScansString();
  precursor_id_ = header_ptr->getPrecId();
  prec_mass_ = header_ptr->getPrecMonoMass();
  spectrum_num_ = spectrum_num;
  proteo_ptr_= proteo_ptr;
  seq_id_ = proteo_ptr->getDbResSeqPtr()->getId();
  seq_name_ = proteo_ptr->getDbResSeqPtr()->getName();
  seq_desc_ = proteo_ptr->getDbResSeqPtr()->getDesc();
  score_ = score;
}

SimplePrsm::SimplePrsm(xercesc::DOMElement* element){
  spectrum_id_ = getIntChildValue(element, "spectrum_id", 0);
  spectrum_scan_ = getChildValue(element, "spectrum_scan", 0);
  precursor_id_ = getIntChildValue(element, "precursor_id", 0);
  prec_mass_ = getDoubleChildValue(element, "precursor_mass", 0);
  spectrum_num_ = getDoubleChildValue(element, "spectrum_number", 0);
  seq_id_ = getIntChildValue(element, "sequence_id", 0);
  seq_name_ = getChildValue(element, "sequence_name", 0);
  seq_desc_ = getChildValue(element, "sequence_desc", 0);
  score_ = getDoubleChildValue(element, "score", 0);
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
  str = convertToString(spectrum_num_);
  xml_doc->addElement(element, "spectrum_number", str.c_str());
  str = convertToString(seq_id_);
  xml_doc->addElement(element, "sequence_id", str.c_str());
  str = convertToString(score_);
  xml_doc->addElement(element, "score", str.c_str());
  xml_doc->addElement(element, "sequence_name", seq_name_.c_str());
  xml_doc->addElement(element, "sequence_desc", seq_desc_.c_str());
  return element;
}

bool SimplePrsm::isSameSpectrum(MsHeaderPtr header_ptr){
  int new_spectrum_id = header_ptr->getId();
  std::string new_spectrum_scan = header_ptr->getScansString();
  int new_precursor_id = header_ptr->getPrecId();
  double new_precursor_mass = header_ptr->getPrecMonoMass();
  if(new_spectrum_id == spectrum_id_ 
     && new_precursor_id == precursor_id_){
    if(std::abs(new_precursor_mass-prec_mass_) > 0.00001 ||
       new_spectrum_scan != spectrum_scan_){
      LOG_ERROR("Error in combine simple Prsms! ");
    }
    return true;
  }
  else {
    return false;
  }
}

bool SimplePrsm::isLargerSpectrumId(MsHeaderPtr header_ptr){
  int new_spectrum_id = header_ptr->getId();
  int new_precursor_id = header_ptr->getPrecId();
  if (new_spectrum_id < spectrum_id_) {
    return true;
  }
  if (new_spectrum_id == spectrum_id_ &&
      new_precursor_id < precursor_id_) {
    return true;
  }
  return false;
}

void SimplePrsm::assignProteoformPtr(const ProteoformPtrVec &proteo_ptrs,
                                     const ProteoformPtrVec2D &mod_proteo_2d_ptrs){
  proteo_ptr_ = proteo_ptrs[seq_id_];
  mod_proteo_ptrs_ = mod_proteo_2d_ptrs[seq_id_];
  if(proteo_ptr_->getSeqId() != seq_id_ || proteo_ptr_->getSeqName() != seq_name_){
    std::cout<< "Sequence ID and/or name is not consistent!" << std::endl;
    std::exit(0);
  }
}

void SimplePrsm::addProteoformPtr(faidx_t *fai, const ResiduePtrVec &residue_list,
                                  const ProtModPtrVec &prot_mods) {
  proteo_ptr_ = readFastaToProteoform(fai, seq_id_, seq_name_, seq_desc_, residue_list); 
  mod_proteo_ptrs_ = generateProtModProteoform(proteo_ptr_, prot_mods);
}


SimplePrsmPtrVec readSimplePrsms(const std::string &file_name){
  SimplePrsmPtrVec result_ptrs;
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if(parser){
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* root = doc.getDocumentElement();
    int simple_prsm_num = prot::getChildCount(root, "simple_prsm");
    for (int i = 0; i < simple_prsm_num; i++) {
      xercesc::DOMElement* simple_prsm = getChildElement(root, "simple_prsm", i);
      result_ptrs.push_back(SimplePrsmPtr(new SimplePrsm(simple_prsm)));
    }
  }
  return result_ptrs;
}

SimplePrsmPtrVec getMatchedSimplePrsms(const SimplePrsmPtrVec &simple_prsm_ptrs,
                                       MsHeaderPtr header_ptr){
  SimplePrsmPtrVec result_ptrs ;
  for(size_t i=0;i<simple_prsm_ptrs.size();i++){
    SimplePrsmPtr prsm_ptr = simple_prsm_ptrs[i];
    if(prsm_ptr->isSameSpectrum(header_ptr)){
      result_ptrs.push_back(prsm_ptr);
    }
  }
  return result_ptrs;
}

SimplePrsmPtrVec getUniqueMatches(SimplePrsmPtrVec &match_ptrs) {
  std::sort(match_ptrs.begin(), match_ptrs.end(),simplePrsmSeqIdUpScoreDown);

  SimplePrsmPtrVec unique_match_ptrs;
  int prev_seq_id = -1;
  for(size_t i=0;i< match_ptrs.size();i++){
    int cur_seq_id = match_ptrs[i]->getSeqId();
    if (cur_seq_id != prev_seq_id) {
      unique_match_ptrs.push_back(match_ptrs[i]);
      prev_seq_id = cur_seq_id;
    }
  }

  return unique_match_ptrs;
}



} /* namespace prot */
