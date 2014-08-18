#include <iostream>
#include <cmath>

#include "simple_prsm.hpp"
#include "base/logger.hpp"
#include "base/proteoform.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

SimplePrsm::SimplePrsm(MsHeaderPtr header_ptr,ProteoformPtr proteo_ptr,int score){
  spectrum_id_ = header_ptr->getId();
  spectrum_scan_ = header_ptr->getScansString();
  precursor_id_ = header_ptr->getPrecId();
  prec_mass_ = header_ptr->getPrecMonoMass();
  proteo_ptr_= proteo_ptr;
  seq_id_ = proteo_ptr->getDbResSeqPtr()->getId();
  seq_name_ = proteo_ptr->getDbResSeqPtr()->getName();
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

} /* namespace prot */
