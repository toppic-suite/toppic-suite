#include "base/logger.hpp"
#include "spec/ms.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

Prsm::Prsm(ProteoformPtr proteoform_ptr, const DeconvMsPtrVec &deconv_ms_ptr_vec, 
           double adjusted_prec_mass, SpParaPtr sp_para_ptr) {
  proteoform_ptr_ = proteoform_ptr;
  deconv_ms_ptr_vec_ = deconv_ms_ptr_vec;
  MsHeaderPtr header_ptr = deconv_ms_ptr_vec[0]->getHeaderPtr();
  spectrum_id_ = header_ptr->getId();
  spectrum_scan_ = header_ptr->getScansString();
  precursor_id_ = header_ptr->getPrecId();
  spectrum_num_ = deconv_ms_ptr_vec.size();
  ori_prec_mass_ = header_ptr->getPrecMonoMass();
  adjusted_prec_mass_ = adjusted_prec_mass;
  //calibration_ = calibration;
  init(sp_para_ptr);
}

void Prsm::init(SpParaPtr sp_para_ptr) {
  refine_ms_three_vec_ = createMsThreePtrVec(deconv_ms_ptr_vec_, sp_para_ptr, adjusted_prec_mass_);
  //refine_ms_three_-> recalibrate(calibration_);
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
  std::sort(pairs.begin(),pairs.end(),peakIonPairUp);
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
  PeakIonPairPtrVec pairs = getPeakIonPairs(proteoform_ptr_, refine_ms_three_vec_,
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
  std::sort(pairs.begin(), pairs.end(), peakIonPairUp);
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
    PeakIonPairPtrVec pairs = getPeakIonPairs (proteoform_ptr_, refine_ms_three_vec_[i], 
                                               sp_para_ptr->getMinMass());
    match_fragment_num_ += compMatchFragNum(pairs);
    match_peak_num_ += compMatchPeakNum(pairs);
  }
}

xercesc::DOMElement* Prsm::toXmlElement(XmlDOMDocument* xml_doc){
	xercesc::DOMElement* element = xml_doc->createElement("prsm");
  //LOG_DEBUG("Element created");
	std::string str = convertToString(prsm_id_);
	xml_doc->addElement(element, "prsm_id", str.c_str());
	str = convertToString(spectrum_id_);
	xml_doc->addElement(element, "spectrum_id", str.c_str());
	xml_doc->addElement(element, "spectrum_scan", spectrum_scan_.c_str());
	str = convertToString(precursor_id_);
	xml_doc->addElement(element, "precursor_id", str.c_str());
	str = convertToString(spectrum_num_);
	xml_doc->addElement(element, "spectrum_number", str.c_str());
	str = convertToString(ori_prec_mass_);
	xml_doc->addElement(element, "ori_prec_mass", str.c_str());
	str = convertToString(adjusted_prec_mass_);
	xml_doc->addElement(element, "adjusted_prec_mass", str.c_str());
	str = convertToString(fdr_);
	xml_doc->addElement(element, "fdr", str.c_str());
	str = convertToString(match_peak_num_);
	xml_doc->addElement(element, "match_peak_num", str.c_str());
	str = convertToString(match_fragment_num_);
	xml_doc->addElement(element, "match_fragment_num", str.c_str());
	proteoform_ptr_->appendXml(xml_doc,element);
	if(prob_ptr_!=nullptr){
	  prob_ptr_->appendXml(xml_doc,element);
	}
  return element;
}

void Prsm::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = toXmlElement(xml_doc);
  parent->appendChild(element);
}

void Prsm::parseXml(xercesc::DOMElement *element) {
  prsm_id_=getIntChildValue(element, "prsm_id", 0);
  spectrum_id_=getIntChildValue(element, "spectrum_id", 0);
  spectrum_scan_=getChildValue(element, "spectrum_scan", 0);
  precursor_id_=getIntChildValue(element, "precursor_id", 0);
  spectrum_num_ = getIntChildValue(element, "spectrum_number", 0);
  ori_prec_mass_=getDoubleChildValue(element, "ori_prec_mass", 0);
  //LOG_DEBUG("reade original precursor mass " << ori_prec_mass_);
  adjusted_prec_mass_=getDoubleChildValue(element, "adjusted_prec_mass", 0);
  fdr_=getDoubleChildValue(element, "fdr", 0);
  match_peak_num_=getDoubleChildValue(element, "match_peak_num", 0);
  match_fragment_num_=getDoubleChildValue(element, "match_fragment_num", 0);

  int prob_count = getChildCount(element,"extreme_value");
  if(prob_count!=0){
    xercesc::DOMElement* prob_element 
        = getChildElement(element,"extreme_value",0);
    prob_ptr_ = ExtremeValuePtr(new ExtremeValue(prob_element));
  }
}

Prsm::Prsm(xercesc::DOMElement* element,ProteoformPtrVec proteoforms){
    parseXml(element);

    xercesc::DOMElement* proteoform_element
        = getChildElement(element,"proteoform",0);
    proteoform_ptr_ 
        = ProteoformPtr(new Proteoform(proteoform_element,proteoforms));
}

Prsm::Prsm(xercesc::DOMElement* element, faidx_t *fai,
           const ResiduePtrVec &residue_ptr_vec) {
  parseXml(element);
  xercesc::DOMElement* proteoform_element
      = getChildElement(element,"proteoform",0);
  proteoform_ptr_ = ProteoformPtr(
      new Proteoform(proteoform_element, fai, residue_ptr_vec));
}

double Prsm::getEValue() {
  if (prob_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  }
  else {
    return prob_ptr_->getEValue();
  }
}

double Prsm::getPValue() {
  if (prob_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  }
  else {
    return prob_ptr_->getPValue();
  }
}

double Prsm::getOneProtProb() {
  if (prob_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  }
  else {
    return prob_ptr_->getOneProtProb();
  }
}


bool Prsm::isMatchMs(MsHeaderPtr header_ptr) {
  int id = header_ptr->getId();
  std::string scan = header_ptr->getScansString();
  int prec_id = header_ptr->getPrecId();
  double prec_mass = header_ptr->getPrecMonoMass();
  if (id == spectrum_id_ && prec_id == precursor_id_) {
    LOG_TRACE("scan " << scan << " prsm scan " << spectrum_scan_
              << " prec mass " << prec_mass << " prsm mass " << ori_prec_mass_);
    if (scan != spectrum_scan_) {
      LOG_ERROR("Error in Prsm.");
    }
    return true;
  } else {
    LOG_TRACE("prsm spectrum id " << spectrum_id_ << " ms spectrum id " << id);
    return false;
  }
}


PrsmPtrVec readPrsm(const std::string &file_name, const ProteoformPtrVec &proteo_ptrs){
    PrsmPtrVec results;
    XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
    if(parser){
        XmlDOMDocument doc(parser, file_name.c_str());
        xercesc::DOMElement* root = doc.getDocumentElement();
        int simple_prsm_num = getChildCount(root, "prsm");
        for (int i = 0; i < simple_prsm_num; i++) {
            xercesc::DOMElement* prsm_element = getChildElement(root, "prsm", i);
            results.push_back(PrsmPtr(new Prsm(prsm_element, proteo_ptrs)));
        }
    }
    return results;
}


void addSpectrumPtrsToPrsms(PrsmPtrVec &prsm_ptrs, PrsmParaPtr prsm_para_ptr){
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(), prsm_para_ptr->getGroupSpecNum());
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  SpectrumSetPtr spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr);
  //use prsm order information (ordered by spectrum id then prec id)
  int start_prsm = 0;
  while (spec_set_ptr != nullptr) {
    if (spec_set_ptr->isValid()) {
      DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
      int spectrum_id = deconv_ms_ptr_vec[0]->getHeaderPtr()->getId();
      int prec_id = deconv_ms_ptr_vec[0]->getHeaderPtr()->getPrecId();
      LOG_TRACE("spectrum id " << spectrum_id);
      for(size_t i = start_prsm; i<prsm_ptrs.size(); i++){
        if(prsm_ptrs[i]->isMatchMs(deconv_ms_ptr_vec[0]->getHeaderPtr())){
          prsm_ptrs[i]->setDeconvMsPtrVec(deconv_ms_ptr_vec);
          double new_prec_mass = prsm_ptrs[i]->getAdjustedPrecMass();
          prsm_ptrs[i]->setRefineMsVec(createMsThreePtrVec(deconv_ms_ptr_vec, prsm_para_ptr->getSpParaPtr(),
                                                           new_prec_mass));
        }
        if ((spectrum_id == prsm_ptrs[i]->getSpectrumId() && prec_id < prsm_ptrs[i]->getPrecursorId()) ||
            spectrum_id < prsm_ptrs[i]->getSpectrumId()) {
          start_prsm = i;
          break;
        }
      }
    }
    spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr);
  }
  reader.close();
}

void filterPrsms(const PrsmPtrVec &prsm_ptrs, MsHeaderPtr header_ptr, 
                 PrsmPtrVec &sele_prsm_ptrs) {
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->isMatchMs(header_ptr)) {
      sele_prsm_ptrs.push_back(prsm_ptrs[i]);
    }
  }
}

std::vector<int> getSpeciesIds(const PrsmPtrVec &prsm_ptrs){
  std::vector<int> species_ids;
  for(size_t i=0;i<prsm_ptrs.size();i++){
    bool find = false;
    for(size_t j=0;j<species_ids.size();j++){
      if(species_ids[j]==prsm_ptrs[i]->getProteoformPtr()->getSpeciesId()){
        find = true;
        break;
      }
    }
    if(!find){
      species_ids.push_back(prsm_ptrs[i]->getProteoformPtr()->getSpeciesId());
    }
  }
  return species_ids;
}

std::vector<int> getSpeciesIds(const PrsmPtrVec &prsm_ptrs,int seq_id){
  std::vector<int> species_ids;
  for(size_t i=0;i<prsm_ptrs.size();i++){
    int new_id = prsm_ptrs[i]->getProteoformPtr()->getSpeciesId();
    if(prsm_ptrs[i]->getProteoformPtr()->getDbResSeqPtr()->getId() == seq_id){
      bool flag= false;
      for(size_t j=0;j<species_ids.size();j++){
        if(species_ids[j]==new_id){
          flag=true;
          break;
        }
      }
      if(!flag){
        species_ids.push_back(new_id);
      }
    }
  }
  return species_ids;
}

PrsmPtrVec selectSpeciesPrsms(const PrsmPtrVec &prsm_ptrs,int species_id){
  PrsmPtrVec select_prsm_ptrs;
  for(size_t i=0;i<prsm_ptrs.size();i++){
    if(species_id == prsm_ptrs[i]->getProteoformPtr()->getSpeciesId()){
      select_prsm_ptrs.push_back(prsm_ptrs[i]);
    }
  }
  return select_prsm_ptrs;
}

}
