#include "base/logger.hpp"
#include "base/proteoform_factory.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "spec/ms.hpp"
#include "spec/extend_ms_factory.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "prsm/peak_ion_pair_factory.hpp"
#include "prsm/prsm.hpp"

namespace prot {

Prsm::Prsm(ProteoformPtr proteoform_ptr, const DeconvMsPtrVec &deconv_ms_ptr_vec, 
           double adjusted_prec_mass, SpParaPtr sp_para_ptr):
    adjusted_prec_mass_(adjusted_prec_mass),
    proteoform_ptr_(proteoform_ptr),
    deconv_ms_ptr_vec_(deconv_ms_ptr_vec) {
      MsHeaderPtr header_ptr = deconv_ms_ptr_vec[0]->getMsHeaderPtr();
      spectrum_id_ = header_ptr->getId();
      spectrum_scan_ = header_ptr->getScansString();
      precursor_id_ = header_ptr->getPrecId();
      spectrum_num_ = deconv_ms_ptr_vec.size();
      ori_prec_mass_ = header_ptr->getPrecMonoMass();
      init(sp_para_ptr);
    }

Prsm::Prsm(xercesc::DOMElement* element, FastaIndexReaderPtr reader_ptr, 
           const ModPtrVec &fix_mod_list) {
  parseXml(element);
  std::string form_elem_name = Proteoform::getXmlElementName();
  xercesc::DOMElement* form_element
      = XmlDomUtil::getChildElement(element, form_elem_name.c_str(),0);
  proteoform_ptr_ = ProteoformPtr(new Proteoform(form_element, reader_ptr, fix_mod_list));
}


void Prsm::init(SpParaPtr sp_para_ptr) {
  refine_ms_three_vec_ 
      = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec_, sp_para_ptr, adjusted_prec_mass_);
  initScores(sp_para_ptr);
}

inline double compMatchFragNum(const PeakIonPairPtrVec &pairs) {
  double match_fragment_num = 0;
  TheoPeakPtr prev_ion(nullptr);;
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i]->getTheoPeakPtr() != prev_ion) {
      prev_ion = pairs[i]->getTheoPeakPtr();
      match_fragment_num += pairs[i]->getRealPeakPtr()->getScore();
    }
  }
  return match_fragment_num;
}

inline double compMatchPeakNum(PeakIonPairPtrVec &pairs) {
  double match_peak_num = 0;
  std::sort(pairs.begin(),pairs.end(),PeakIonPair::cmpRealPeakPosInc);
  DeconvPeakPtr prev_deconv_peak(nullptr);
  //  LOG_DEBUG("total peak number " << pairs.size());
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i]->getRealPeakPtr()->getBasePeakPtr() != prev_deconv_peak) {
      prev_deconv_peak = pairs[i]->getRealPeakPtr()->getBasePeakPtr();
      match_peak_num += pairs[i]->getRealPeakPtr()->getScore();
    }
  }
  return match_peak_num;
}

void Prsm::initMatchNum(double min_mass) {
  PeakIonPairPtrVec pairs = 
      PeakIonPairFactory::genePeakIonPairs(proteoform_ptr_, refine_ms_three_vec_,
                                           min_mass);
  match_peak_num_ = 0;
  match_fragment_num_ = 0;
  TheoPeakPtr prev_ion(nullptr);
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i]->getTheoPeakPtr() != prev_ion) {
      prev_ion = pairs[i]->getTheoPeakPtr();
      match_fragment_num_ += pairs[i]->getRealPeakPtr()->getScore();
    }
  }
  std::sort(pairs.begin(), pairs.end(), PeakIonPair::cmpRealPeakPosInc);
  DeconvPeakPtr prev_deconv_peak(nullptr);
  for (size_t i = 0; i < pairs.size(); i++) {
    if (pairs[i]->getRealPeakPtr()->getBasePeakPtr() != prev_deconv_peak) {
      prev_deconv_peak = pairs[i]->getRealPeakPtr()->getBasePeakPtr();
      match_peak_num_ += pairs[i]->getRealPeakPtr()->getScore();
    }
  }
}

void Prsm::initScores(SpParaPtr sp_para_ptr) {
  match_fragment_num_ = 0;
  match_peak_num_ = 0;
  for (size_t i = 0; i < refine_ms_three_vec_.size(); i++) {
    // refined one 
    PeakIonPairPtrVec pairs = 
        PeakIonPairFactory::genePeakIonPairs(proteoform_ptr_, refine_ms_three_vec_[i], 
                                             sp_para_ptr->getMinMass());
    match_fragment_num_ += compMatchFragNum(pairs);
    match_peak_num_ += compMatchPeakNum(pairs);
  }
}

xercesc::DOMElement* Prsm::toXmlElement(XmlDOMDocument* xml_doc){
  std::string element_name = Prsm::getXmlElementName();
	xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
	std::string str = StringUtil::convertToString(prsm_id_);
	xml_doc->addElement(element, "prsm_id", str.c_str());
	str = StringUtil::convertToString(spectrum_id_);
	xml_doc->addElement(element, "spectrum_id", str.c_str());
	xml_doc->addElement(element, "spectrum_scan", spectrum_scan_.c_str());
	str = StringUtil::convertToString(precursor_id_);
	xml_doc->addElement(element, "precursor_id", str.c_str());
	str = StringUtil::convertToString(spectrum_num_);
	xml_doc->addElement(element, "spectrum_number", str.c_str());
	str = StringUtil::convertToString(ori_prec_mass_);
	xml_doc->addElement(element, "ori_prec_mass", str.c_str());
	str = StringUtil::convertToString(adjusted_prec_mass_);
	xml_doc->addElement(element, "adjusted_prec_mass", str.c_str());
	str = StringUtil::convertToString(fdr_);
	xml_doc->addElement(element, "fdr", str.c_str());
	str = StringUtil::convertToString(match_peak_num_);
	xml_doc->addElement(element, "match_peak_num", str.c_str());
	str = StringUtil::convertToString(match_fragment_num_);
	xml_doc->addElement(element, "match_fragment_num", str.c_str());
	proteoform_ptr_->appendXml(xml_doc,element);
	if(extreme_value_ptr_!=nullptr){
	  extreme_value_ptr_->appendXml(xml_doc,element);
	}
  return element;
}

void Prsm::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = toXmlElement(xml_doc);
  parent->appendChild(element);
}

void Prsm::parseXml(xercesc::DOMElement *element) {
  prsm_id_=XmlDomUtil::XmlDomUtil::getIntChildValue(element, "prsm_id", 0);
  spectrum_id_=XmlDomUtil::getIntChildValue(element, "spectrum_id", 0);
  spectrum_scan_=XmlDomUtil::getChildValue(element, "spectrum_scan", 0);
  precursor_id_=XmlDomUtil::getIntChildValue(element, "precursor_id", 0);
  spectrum_num_ = XmlDomUtil::getIntChildValue(element, "spectrum_number", 0);
  ori_prec_mass_=XmlDomUtil::getDoubleChildValue(element, "ori_prec_mass", 0);
  adjusted_prec_mass_=XmlDomUtil::getDoubleChildValue(element, "adjusted_prec_mass", 0);
  fdr_=XmlDomUtil::getDoubleChildValue(element, "fdr", 0);
  match_peak_num_=XmlDomUtil::getDoubleChildValue(element, "match_peak_num", 0);
  match_fragment_num_=XmlDomUtil::getDoubleChildValue(element, "match_fragment_num", 0);

  int prob_count = XmlDomUtil::getChildCount(element,"extreme_value");
  if(prob_count!=0){
    xercesc::DOMElement* prob_element 
        = XmlDomUtil::getChildElement(element,"extreme_value",0);
    extreme_value_ptr_ = ExtremeValuePtr(new ExtremeValue(prob_element));
  }
}

double Prsm::getEValue() {
  if (extreme_value_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  }
  else {
    return extreme_value_ptr_->getEValue();
  }
}

double Prsm::getPValue() {
  if (extreme_value_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  }
  else {
    return extreme_value_ptr_->getPValue();
  }
}

double Prsm::getOneProtProb() {
  if (extreme_value_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  }
  else {
    return extreme_value_ptr_->getOneProtProb();
  }
}


// sort by number of matched fragment ions, then start position 
bool Prsm::cmpMatchFragDecStartPosInc(const PrsmPtr &a, const PrsmPtr &b) {
  if(a->getMatchFragNum() > b->getMatchFragNum()){
    return true;
  }
  else if(a->getMatchFragNum() == b->getMatchFragNum()){
    return a->getProteoformPtr()->getStartPos()<b->getProteoformPtr()->getStartPos();
  }
  return false;
}

// sort by the order of spectrum id, the precursor id
bool PrsmCmpSpectrumIdIncPrecursorIdInc(const PrsmPtr &a, const PrsmPtr &b){
  if(a->getSpectrumId() < b->getSpectrumId()){
    return true;
  }
  else if(a->getSpectrumId() > b->getSpectrumId()){
    return false;
  }
  else{
    if(a->getPrecursorId() < b->getPrecursorId()){
      return true;
    }
    return false;
  }
}

// sort by spectrum id then match ions
bool Prsm::cmpSpectrumIdIncMatchFragDec(const PrsmPtr &a, const PrsmPtr &b){
  if(a->getSpectrumId() < b->getSpectrumId()){
    return true;
  }
  else if(a->getSpectrumId() > b->getSpectrumId()){
    return false;
  }
  else{
    if(a->getMatchFragNum() > b->getMatchFragNum()){
      return true;
    }
    return false;
  }
}

// sort by spectrum id then evalue
bool Prsm::cmpSpectrumIdIncEvalueInc(const PrsmPtr &a, const PrsmPtr &b){
  if(a->getSpectrumId() < b->getSpectrumId()){
    return true;
  }
  else if(a->getSpectrumId() > b->getSpectrumId()){
    return false;
  }
  else{
    if(a->getEValue() < b->getEValue()){
      return true;
    }
    return false;
  }
}

}
