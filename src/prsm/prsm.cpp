//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "ms/spec/extend_ms_factory.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm.hpp"

namespace toppic {

Prsm::Prsm(ProteoformPtr proteoform_ptr, const DeconvMsPtrVec &deconv_ms_ptr_vec,
           double adjusted_prec_mass, SpParaPtr sp_para_ptr):
    adjusted_prec_mass_(adjusted_prec_mass),
    proteoform_ptr_(proteoform_ptr),
    deconv_ms_ptr_vec_(deconv_ms_ptr_vec) {
      MsHeaderPtr header_ptr = deconv_ms_ptr_vec[0]->getMsHeaderPtr();
      spectrum_id_ = header_ptr->getId();
      spectrum_scan_ = header_ptr->getScansString();
      precursor_id_ = header_ptr->getPrecId();
//    prec_feature_id_ = header_ptr->getFeatureId();
//    prec_feature_inte_ = header_ptr->getFeatureInte();
      spectrum_num_ = deconv_ms_ptr_vec.size();
      ori_prec_mass_ = header_ptr->getPrecMonoMass();
      init(sp_para_ptr);
    }

Prsm::Prsm(XmlDOMElement* element, FastaIndexReaderPtr reader_ptr,
           const ModPtrVec &fix_mod_list) {
  parseXml(element);
  std::string form_elem_name = Proteoform::getXmlElementName();
  XmlDOMElement* form_element
      = xml_dom_util::getChildElement(element, form_elem_name.c_str(), 0);
  proteoform_ptr_ = std::make_shared<Proteoform>(form_element, reader_ptr, fix_mod_list);
}

Prsm::Prsm(const Prsm &obj) {
  adjusted_prec_mass_ = obj.adjusted_prec_mass_;
  proteoform_ptr_ = obj.proteoform_ptr_;
  deconv_ms_ptr_vec_ = obj.deconv_ms_ptr_vec_;
  file_name_ = obj.file_name_;
  spectrum_id_ = obj.spectrum_id_;
  spectrum_scan_ = obj.spectrum_scan_;
  precursor_id_ = obj.precursor_id_;
  prec_feature_id_ = obj.prec_feature_id_;
  prec_feature_inte_ = obj.prec_feature_inte_;
  spectrum_num_ = obj.spectrum_num_;
  ori_prec_mass_ = obj.ori_prec_mass_;
  match_peak_num_ = obj.match_peak_num_;
  match_fragment_num_ = obj.match_fragment_num_;
}

void Prsm::init(SpParaPtr sp_para_ptr) {
  refine_ms_three_vec_
      = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec_, sp_para_ptr, adjusted_prec_mass_);
  initScores(sp_para_ptr);
}

void Prsm::initMatchNum(double min_mass) {
  PeakIonPairPtrVec pairs =
      peak_ion_pair_util::genePeakIonPairs(proteoform_ptr_, refine_ms_three_vec_, min_mass);

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
        peak_ion_pair_util::genePeakIonPairs(proteoform_ptr_, refine_ms_three_vec_[i],
                                             sp_para_ptr->getMinMass());
    match_fragment_num_ += peak_ion_pair_util::compMatchFragNum(pairs);
    match_peak_num_ += peak_ion_pair_util::compMatchPeakNum(pairs);
  }
}

XmlDOMElement* Prsm::toXmlElement(XmlDOMDocument* xml_doc) {
  std::string element_name = Prsm::getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "file_name", file_name_.c_str());
  std::string str = str_util::toString(prsm_id_);
  xml_doc->addElement(element, "prsm_id", str.c_str());
  str = str_util::toString(spectrum_id_);
  xml_doc->addElement(element, "spectrum_id", str.c_str());
  xml_doc->addElement(element, "spectrum_scan", spectrum_scan_.c_str());
  str = str_util::toString(precursor_id_);
  xml_doc->addElement(element, "precursor_id", str.c_str());
  str = str_util::toString(prec_feature_id_);
  xml_doc->addElement(element, "precursor_feature_id", str.c_str());
  str = str_util::toString(prec_feature_inte_);
  xml_doc->addElement(element, "precursor_feature_inte", str.c_str());
  str = str_util::toString(frac_feature_score_);
  xml_doc->addElement(element, "frac_feature_score", str.c_str());
  str = str_util::toString(spectrum_num_);
  xml_doc->addElement(element, "spectrum_number", str.c_str());
  str = str_util::toString(ori_prec_mass_);
  xml_doc->addElement(element, "ori_prec_mass", str.c_str());
  str = str_util::toString(adjusted_prec_mass_);
  xml_doc->addElement(element, "adjusted_prec_mass", str.c_str());
  str = str_util::toString(fdr_);
  xml_doc->addElement(element, "fdr", str.c_str());
  str = str_util::toString(proteoform_fdr_);
  xml_doc->addElement(element, "proteoform_fdr", str.c_str());
  str = str_util::toString(match_peak_num_);
  xml_doc->addElement(element, "match_peak_num", str.c_str());
  str = str_util::toString(match_fragment_num_);
  xml_doc->addElement(element, "match_fragment_num", str.c_str());
  str = str_util::toString(getNormMatchFragNum());
  xml_doc->addElement(element, "norm_match_fragment_num", str.c_str());
  proteoform_ptr_->appendXml(xml_doc, element);
  if (extreme_value_ptr_ != nullptr) {
    extreme_value_ptr_->appendXml(xml_doc, element);
  }
  return element;
}

void Prsm::appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  XmlDOMElement* element = toXmlElement(xml_doc);
  parent->appendChild(element);
}

void Prsm::parseXml(XmlDOMElement *element) {
  file_name_ = xml_dom_util::getChildValue(element, "file_name", 0);
  prsm_id_ = xml_dom_util::getIntChildValue(element, "prsm_id", 0);
  spectrum_id_ = xml_dom_util::getIntChildValue(element, "ms/spectrum_id", 0);
  spectrum_scan_ = xml_dom_util::getChildValue(element, "ms/spectrum_scan", 0);
  precursor_id_ = xml_dom_util::getIntChildValue(element, "precursor_id", 0);
  prec_feature_id_ = xml_dom_util::getIntChildValue(element, "precursor_feature_id", 0);
  prec_feature_inte_ = xml_dom_util::getDoubleChildValue(element, "precursor_feature_inte", 0);
  frac_feature_score_ = xml_dom_util::getDoubleChildValue(element, "frac_feature_score", 0);
  spectrum_num_ = xml_dom_util::getIntChildValue(element, "ms/spectrum_number", 0);
  ori_prec_mass_ = xml_dom_util::getDoubleChildValue(element, "ori_prec_mass", 0);
  adjusted_prec_mass_ = xml_dom_util::getDoubleChildValue(element, "adjusted_prec_mass", 0);
  fdr_ = xml_dom_util::getDoubleChildValue(element, "fdr", 0);
  proteoform_fdr_ = xml_dom_util::getDoubleChildValue(element, "proteoform_fdr", 0);
  match_peak_num_ = xml_dom_util::getDoubleChildValue(element, "match_peak_num", 0);
  match_fragment_num_ = xml_dom_util::getDoubleChildValue(element, "match_fragment_num", 0);

  int prob_count = xml_dom_util::getChildCount(element, "extreme_value");
  if (prob_count != 0) {
    XmlDOMElement* prob_element
        = xml_dom_util::getChildElement(element, "extreme_value", 0);
    extreme_value_ptr_ = std::make_shared<ExtremeValue>(prob_element);
  }
}

double Prsm::getEValue() {
  if (extreme_value_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  } else {
    return extreme_value_ptr_->getEValue();
  }
}

double Prsm::getPValue() {
  if (extreme_value_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  } else {
    return extreme_value_ptr_->getPValue();
  }
}

double Prsm::getOneProtProb() {
  if (extreme_value_ptr_ == nullptr) {
    LOG_WARN("Probability pointer is null.");
    return -1;
  } else {
    return extreme_value_ptr_->getOneProtProb();
  }
}

/* this function is tempory for testing mass graph alignment */
double Prsm::getNormMatchFragNum() {
  int var_change_num = proteoform_ptr_->getVariablePtmNum();

  int unexp_change_num = proteoform_ptr_->getMassShiftNum(AlterType::UNEXPECTED);

  int start_pos = proteoform_ptr_->getStartPos();

  int end_pos = proteoform_ptr_->getEndPos();

  double score = match_fragment_num_ - 2 * var_change_num - 6 * unexp_change_num;

  int trunc_len = getProteoformPtr()->getProtModPtr()->getTruncPtr()->getTruncLen();

  if (start_pos == trunc_len) {
    score += 1;
  }

  if (end_pos == getProteoformPtr()->getFastaSeqPtr()->getAcidPtmPairLen() - 1) {
    score += 1;
  }

  if (score < 0) score = 0.0;

  return score;
}

// sort by the number of matched fragments, then the number of matched peaks
bool Prsm::cmpMatchFragmentDecMatchPeakDec(const PrsmPtr &a, const PrsmPtr &b) {
  if (a->getMatchFragNum() > b->getMatchFragNum()) {
    return true;
  } else if (a->getMatchFragNum() < b->getMatchFragNum()) {
    return false;
  }
  return a->getMatchPeakNum() > b->getMatchPeakNum();
}

// sort by number of matched fragment ions, then start position
bool Prsm::cmpMatchFragDecStartPosInc(const PrsmPtr &a, const PrsmPtr &b) {
  if (a->getMatchFragNum() > b->getMatchFragNum()) {
    return true;
  } else if (a->getMatchFragNum() == b->getMatchFragNum()) {
    return a->getProteoformPtr()->getStartPos() < b->getProteoformPtr()->getStartPos();
  }
  return false;
}

// sort by the order of spectrum id, the precursor id
bool Prsm::cmpSpectrumIdIncPrecursorIdInc(const PrsmPtr &a, const PrsmPtr &b) {
  if (a->getSpectrumId() < b->getSpectrumId()) {
    return true;
  } else if (a->getSpectrumId() > b->getSpectrumId()) {
    return false;
  } else {
    if (a->getPrecursorId() < b->getPrecursorId()) {
      return true;
    }
    return false;
  }
}

// sort by spectrum id then match ions
bool Prsm::cmpSpectrumIdIncMatchFragDec(const PrsmPtr &a, const PrsmPtr &b) {
  if (a->getSpectrumId() < b->getSpectrumId()) {
    return true;
  } else if (a->getSpectrumId() > b->getSpectrumId()) {
    return false;
  } else {
    if (a->getMatchFragNum() > b->getMatchFragNum()) {
      return true;
    }
    return false;
  }
}

// sort by spectrum id then evalue
bool Prsm::cmpSpectrumIdIncEvalueInc(const PrsmPtr &a, const PrsmPtr &b) {
  if (a->getSpectrumId() < b->getSpectrumId()) {
    return true;
  } else if (a->getSpectrumId() > b->getSpectrumId()) {
    return false;
  } else if (a->getEValue() < b->getEValue()) {
    return true;
  } else if (a->getEValue() > b->getEValue()) {
    return false;
  } else {
    return a->getProteoformPtr()->getSeqName() < b->getProteoformPtr()->getSeqName();
  }
}

bool Prsm::cmpNormMatchFragmentDec(const PrsmPtr &a, const PrsmPtr &b) {
  if (a->getNormMatchFragNum() == b->getNormMatchFragNum()) {
    return a->getProteoformPtr()->getVariablePtmNum() < b->getProteoformPtr()->getVariablePtmNum();
  } else {
    return a->getNormMatchFragNum() > b->getNormMatchFragNum();
  }
}

void Prsm::setProteoformPtr(ProteoformPtr proteoform, SpParaPtr sp_para_ptr) {
  proteoform_ptr_ = proteoform;
  init(sp_para_ptr);
}

}  // namespace toppic
